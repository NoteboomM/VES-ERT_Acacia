# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:35:30 2024

@author: matthew.noteboom

This is a script to import, check and invert resistivity data, currently mostly 
a work in progress. Initial testing based on ERT profiles recorded with the 
Supersting instrument in Lubango, Angola (proj 201025). Although Thomas Gunther
encouraged loading the stg files directly, I prefer converting to the cleaner 
'.ohm' pybert format, and importing those (see stg-to-ohm.py). Importing directly
from stg file means we need to separately import the terrain (.trn file).

24/4/24: conversion tool and this script can store/import total IP, but not the 
decay. Have not tested inversion. Inversion of resistivity is stable, but at 
least with the Lubango data requires some QC on spurious measurements (see below).

Also tested with .ohm file converted from Acacia Terrameter .txt file.

August 2024 - testing with layer constraints and 'prior' data

"""
# import useful/necessary packages
import os, sys

import pygimli as pg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pygimli.physics import ert
import pygimli.meshtools as mt
from pygimli.frameworks import PriorModelling, JointModelling
from pygimli.viewer.mpl import draw1DColumn
# import pybert as pb

pyversion = sys.version

# This version works with a temporary testing folder, but can be changed to any path
# (although eventually there is probably a length limit).
os.chdir("C:\\Temp\\211269_OptHof_MSNTesting\\2024-03-07 OPTHOF-ROLL_LS214100268_19-21-11")

#%%%
"""

The '.ohm' format used in one of the pyGIMLi examples seems simple enough: a section
with x and z coords for electrodes, and a section with electrode locations as index
numbers from the list at the top of the file, plus columns for i,u,r,ip,err etc.

The Res2dInv dat format is mostly simple, but I don't really know what's in the
header block.

Ideally, IP processing would include some full-decay characterisation, but in the 
absence of that, a decay check would be good for QC.

"""

data = ert.load("2024-03-07 OPTHOF-ROLL_GradientXL_2.ohm") # Load data to container from BERT format
# data = ert.load("2022-11-15 ZTT BRON21 T0_Gradient_697i_1933p_Acacia_1.ohm") # Load data to container from BERT format
terrain = pg.z(data) # store the terrain information which gets moved to 'y' axis during inversion

# load the table of estimated/measured layers & properties
prior = pd.read_csv("layers.txt", sep = "\\s+")#, header = 0)

print(data, data.tokenList()) # print summary - good to check # of sensors, measurements and channels
# plot terrain profile. set_aspect method defines vertical exaggeration
# fig, ax = plt.subplots()
# ax.plot(pg.x(data), pg.y(data), '.-')
# ax.set_aspect(2.0)

#%%
# calculate geometric factors: False=the simple way using formula with a,b,m,n,
# otherwise 'True' calculates geometric factor numerically from terrain information.
# difference minor in smooth/flat terrain

# Note: Terrameter examples have no terrain(?)
data['k'] = ert.createGeometricFactors(data, numerical=False)
k0 = ert.createGeometricFactors(data, numerical=True) 

# plot geometry factor (multiplied by -1)
ert.show(data, data['k'], logScale=True, label='Geometry factor')

# make a plot of the difference between two approaches
ert.show(data, vals=100*k0/data['k'], label='Topography effect %',
        cMap="bwr", logScale=False); #cMin=0.8, cMax=1.25, 

# calculate resistivity from resistance (v/i) channel and calculated geometric factor, then plot
# but for Terrameter LS file, rhoa already calculated, and with no terrain, no new
# k necessary
data['rhoa_calc'] = data['r'] * data['k']

# want to start labelling and saving some of these plots...but I can't make it work so far!
# fig = pg.plt.figure()
# fig.suptitle('Observed apparent resistivity') # Set title
ert.show(data, data['rhoa'])
ert.show(data, data['rhoa_calc'])
ig = plt.gcf()
fig.suptitle("Observed Apparent Resistivity")
fig.savefig('OptHof_Roll_ObsAppRes.png', dpi=200)

plotkw = dict(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r", cMin=0.2, cMax=200, logScale=True)

invkw = dict(data=data, verbose=True, lam=2, zWeight=0.5, robust=True) # lambda is a smoothness constraint, lower zweight emphasises layering

#%% 
"""
********************
  OPTIONAL BLOCKS!
********************

# these are a couple of filter statements I've found the Lubango data sometimes 
# needs, to eliminate spurious very low (or negative) or very high values. But 
# hard to generalize as true resistivity can vary significantly between locations.
data.remove(data["rhoa"] <= 1500)
data.remove(data["rhoa"] > 800000)

# Also a routine for checking points with nearest neighbours. Except for the last
# value, compares value n with average of n-1 and n+1 (for first value, n-1 takes
# last value)
data['check_adj'] = np.zeros(len(data['r']))
for (n,rhoa) in enumerate(data['rhoa']):
    if n==len(data['rhoa'])-1:
        data['check_adj'][n] = data['rhoa'][n]/data['rhoa'][n-1]
    else:
        data['check_adj'][n] = data['rhoa'][n]/np.mean([data['rhoa'][n-1],data['rhoa'][n+1]])

# plot check field
ert.show(data, data['check_adj'], label='data check field: n/(mean of neighbours)')

# plot the check field, mostly for interest/curiosity
data.remove(data['check_adj'] > 10)
data.remove(data['check_adj'] < 0.01)

# plot resistivity data again after filtering
ert.show(data)
"""
#%%

# data['err'] = ert.estimateError(data, # simulate an error field if necessary
#                                 absoluteUError=0.00005, # 50µV
#                                 relativeError=0.03)  # 3% 

# Set a minimum error of 0.5% (tiny values mess up the misfit calc/plot later)
# BUT!!!! using this seems to confuse the ERTManager in the next block;
# Remove until review
# data["err"][data["err"] < 0.005] = 0.005

# plot a pseudosection of the error
ert.show(data, data['err']*100, label="error [%]", logScale=False)

# Option to remove points for high noise
# data.remove(data["err"] >= 5/100 ) # filter extremes?

#%%
# call the ERTManager class
mgr = ert.ERTManager(data)
# run the inversion, with some settings for mesh, lambda, smoothness
mgr.invert(quality=34, size=1., **invkw)

#%%

# plot comparison pseudosections of data & response of inverted model. If inversion
# worked well, they're indistinguishable. Misfit more useful...
mgr.showFit()
# fig = plt.gcf()
# fig.suptitle("Comparison of data & inversion response")
# fig.savefig('OptHof_ObsAppRes.png', dpi=200)

data['misfit'] = pg.log(mgr.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, data['misfit'], cMap="bwr", cMin = -100, cMax = 100,
        label='misfit: response/data/err',
        ) #, cMin=-10, cMax=10

# plot inverted section, can change cMap from default 'Spectral_r' (maybe test 
# viridis or similar for accessibility), set colour max/min, probably specify log 
# scale, and toy with coverage?
ax, _ = mgr.showResult(**plotkw)
                       # label='Apparent Resistivity ohm.m, zwt=0.1, lambda=2, NOT robust')#, coverage=.5) # , cMap="viridis") #, cMap="Spectral_r"); #cMin=1, cMax=500, 
ax.set_ylim(-25,0)
ax.set_xlim(0,140)
ax.grid()
fig = plt.gcf()
fig.suptitle("Inverted Resistivity model, zwt=0.5, lambda=2, NOT robust", y=0.8)
fig.savefig('Unconstrained_Inversion.png', dpi=200)

#%%

# A little tester here using a manually created mesh with lines to encourage a 
# layered earth model. Create a 'world' that's more extensive than the area
# of interest/surveyed, then add two lines. In this case starting x=25 due to 
# feature near the beginning of the line that doesn't look like layered earth.

world = mt.createWorld(start=[-50, 0.], end=[190, -50], worldMarker=True)
# pg.show(world)

line1 = mt.createLine(start=[25, -1], end=[140, -1], marker=2) # marker>0 means it functions as a constraint

line2 = mt.createLine(start=[25, -10], end=[140, -10], marker=3) # marker>0 means it functions as a constraint

world += line1 + line2

pg.show(world)

#%%

# Load up the list of electrode locations and add nodes to the 'world'...
eleclist = list(set(pg.x(data)))
electrodes = np.zeros((len(pg.x(data)),3), dtype=float)
for n,e in enumerate(eleclist): electrodes[n] = [e, 0., 0.]

for p in electrodes: 
    world.createNode(p)
    world.createNode(p -[0, 0.1, 0])

# ...then create and display a mesh from the world with lines and nodes.
mesh = pg.meshtools.createMesh(world, quality=34, size=1.)

ax, _ = pg.show(mesh)
ax.set_ylim(-25,0)
ax.set_xlim(0,140)
fig = plt.gcf()
fig.suptitle("Custom world & mesh for Op 't Hof inversion", y=0.7)
fig.savefig('LayerMesh.png', dpi=200)

#%%
# Run the inversion on the 'layered' mesh
mgrlayers = ert.ERTManager()
mgrlayers.invert(mesh=mesh, **invkw)

#%%

# View results of layered inversion
mgrlayers.showFit()

ax, _ = mgrlayers.showResult(**plotkw)
                       # label='Apparent Resistivity ohm.m, zwt=0.2, lambda=10, robust');
ax.set_ylim(-25,0)
ax.set_xlim(0,140)
ax.grid()
fig = plt.gcf()
fig.suptitle("Inverted Resistivity model, zwt=0.5, lambda=2, NOT robust", y=0.8)
fig.savefig('Inversion_with_layer_boundaries.png', dpi=200)

misfit_layers = pg.log(mgrlayers.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, misfit_layers, cMap="bwr", cMin=-100, cMax=100, label='misfit: response/data/err')

#%%
"""
Now...let's try using the "prior" data...
This is adapted, with very little change, from the plot_5_ert_with_priors 
example from the pyGimli team.

"""
# replaced by plotkw earlier
# kw = dict(cMin=0.2, cMax=200, logScale=True, cMap="Spectral_r", 
#           xlabel="x (m)", ylabel="z (m)")

# Split up layer data table
x,z,r = prior['x'], prior['z'], prior['r']

# display earlier unconstrained inversion result with column of layer data
ax, cb = mgr.showResult(**plotkw)
zz = np.abs(z)
iz = np.argsort(z)
dz = np.diff(zz[iz])
thk = np.hstack([dz, dz[-1]])
ztop = -zz[iz[0]] - dz[0]/2
colkw = dict(x=x[0], val=r[iz], thk=thk, width=3, ztopo=ztop)
draw1DColumn(ax, **colkw, **plotkw)
ax.grid(True)
ax.set_ylim(-25,0)
ax.set_xlim(0,140)
fig = plt.gcf()
fig.suptitle("Inverted Resistivity model, unconstrained, with a priori column", y=0.8)
fig.savefig('Unconstrained_Inversion+apriori.png', dpi=200)


# We want to extract the resistivity from the mesh at the points where 
# the prior data are available. To this end, we create a list of points
# (pg.Pos class) and use a forward operator that picks the values from the
# model vector according to the cell where each point is located. (See the 
# regularization tutorial for details about that.)
posVec = [pg.Pos(pos) for pos in zip(x, z)]
para = pg.Mesh(mgr.paraDomain)  # make a copy
para.setCellMarkers(pg.IVector(para.cellCount()))
fopDP = PriorModelling(para, posVec)

# Now use that forward operator to extract, store and plot the model values
# compared to the a priori values
fig, ax = plt.subplots()
ax.semilogx(r, z, label="a priori")
resSmooth = fopDP(mgr.model)
ax.semilogx(resSmooth, z, label="ERT smooth")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid(True)
ax.legend()
fig = plt.gcf()
fig.suptitle("Res-depth profile, unconstrained inversion and a priori data")
fig.savefig('Res-depth_ERTSmooth+apriori.png', dpi=200)

# "As alternative to smoothness, we can use a geostatistic model. The 
# vertical range can be well estimated from the DP data using a variogram 
# analysis, we guess 8m. For the horizontal one, we can only guess a ten 
# times higher value.
mgr.inv.setRegularization(2, correlationLengths=[40, 4])
mgr.invert()

ax, cb = mgr.showResult(**plotkw)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)
draw1DColumn(ax, **colkw, **plotkw)
resGeo = fopDP(mgr.model)

# Plot compare 'ground truth' with different inversions so far...
fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
ax.semilogx(resSmooth, z, label="ERT smooth")
#ax.semilogx(res2, z, label="ERT aniso")
ax.semilogx(resGeo, z, label="ERT geostat")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid()
ax.legend()

# Only shows subtle changes. An alternative is using the interfaces as
# constraints as shown earlier. However, we really want to use the 'ground
# truth' data as inputs for inversion. This is easily accomplished by taking 
# the mapping operator that we already use for interpolation as a forward 
# operator.
# 
# We set up an inversion with this mesh, logarithmic transformations and 
# invert the model.
inv = pg.Inversion(fop=fopDP, verbose=True)
inv.mesh = para
tLog = pg.trans.TransLog()
inv.modelTrans = tLog
inv.dataTrans = tLog
inv.setRegularization(correlationLengths=[40, 4])
model = inv.run(r, relativeError=0.2)
ax, cb = pg.show(para, model, **plotkw)
ax.set_ylim(-25,0)
ax.set_xlim(0,140)
draw1DColumn(ax, **colkw, **plotkw)

# We now use the framework JointModelling to combine the ERT and the DP 
# forward operators. So we set up a new ERT modelling operator and join it 
# with fopDP.
fopJoint = JointModelling([mgr.fop, fopDP])
# fopJoint.setMesh(para)
fopJoint.setData([data, pg.Vector(r)])  # needs to have .size() attribute! (?)

# We first test the joint forward operator. We create a modelling vector of 
# constant resistivity and distribute the model response into the two parts 
# that can be looked at individually.
model = pg.Vector(para.cellCount(), 100)
response = fopJoint(model)
respERT = response[:data.size()]
respDP = response[data.size():]
print(respDP)

# The Jacobian can be created and looked up by (wish I knew what it meant!)
fopJoint.createJacobian(model)  # works
J = fopJoint.jacobian()
print(type(J))  # wrong type
ax, cb = pg.show(J)
print(J.mat(0))
ax, cb = pg.show(J.mat(1), markersize=4)

#%%%

# For the joint inversion, concatenate the data and error vectors, create 
# a new inversion instance, set logarithmic transformations and run the 
# inversion.
dataVec = np.concatenate((data["rhoa"], r))
errorVec = np.concatenate((data["err"], np.ones_like(r)*0.2))
inv = pg.Inversion(fop=fopJoint, verbose=True)
transLog = pg.trans.TransLog()
inv.modelTrans = transLog
inv.dataTrans = transLog
inv.run(dataVec, errorVec, startModel=model)
ax, cb = pg.show(para, inv.model, **plotkw)
draw1DColumn(ax, **colkw, **plotkw)
ax.set_ylim(-25,0)
ax.set_xlim(0,140)

# MAYBE there's some improvement? Hard to be sure. Next use geostatistics to
# extend the influence of the hole data.
inv.setRegularization(2, correlationLengths=[40, 4])
model = inv.run(dataVec, errorVec, startModel=model)
ax, cb = pg.show(para, model, **plotkw)
draw1DColumn(ax, **colkw, **plotkw)
ax.set_ylim(-25,0)
ax.set_xlim(0,140)

# We split the model response in the ERT part and the prior data part.
# The first is shown as misfit.
respERT = inv.response[:data.size()]
misfit = - respERT / data["rhoa"] * 100 + 100
ax, cb = ert.show(data, misfit, cMap="bwr", cMin=-5, cMax=5)

# The second is shown as a depth profile.
respDP = inv.response[data.size():]
fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
ax.semilogx(respDP, z, label="response")
ax.semilogx(resSmooth, z, label="ERT smooth")
ax.semilogx(resGeo, z, label="ERT geostat")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid()
ax.legend()


"""

fig.suptitle(sitecode) # Set title
fig.savefig(outname+'.png', dpi=200)


eleclist = list(set(pg.x(data)))
electrodes = np.zeros((len(pg.x(data)),3), dtype=float)
for n,e in enumerate(eleclist):
    electrodes[n] = [e, 0., 0.]

# create 'world' in which we'll create a mesh and run inversion
spacing = 2.

world = mt.createWorld(start=[min(electrodes)-5*spacing, 0.], 
                       end=[max(electrodes)+5*spacing, -10*spacing],
                       worldMarker=True)
pg.show(world)

line1 = mt.createLine(start=[min(electrodes)-4*spacing, -1],
                                end=[max(electrodes)+4*spacing, -1],
                                marker=2) # marker>0 means it functions as a constraint
line2 = mt.createLine(start=[min(electrodes)-4*spacing, -13],
                                end=[max(electrodes)+4*spacing, -13],
                                marker=3) # marker>0 means it functions as a constraint

world += line1 + line2

pg.show(world)

for p in electrodes:
    world.createNode(p)
    # world.createNode(p - [0, 0.5])

invmesh = mt.createMesh(world, quality=34)
ax, _ = pg.show(invmesh)
"""