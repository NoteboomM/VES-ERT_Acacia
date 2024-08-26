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

prior = pd.read_csv("layers.txt", sep = "\\s+")#, header = 0)
# prior['x'] = 100

"""
Alternatively, with thanks to Thomas Gunther, use pybert to import directly from STG file,
but then we need to load the trn file separately (which is not coded as of 27/3/2024):
import pybert as pb
data = pb.importData("???.stg")
print(data, data.tokenList())

"""

print(data, data.tokenList()) # print summary - good to check # of sensors, measurements and channels
# plot terrain profile. set_aspect method defines vertical exaggeration
fig, ax = plt.subplots()
ax.plot(pg.x(data), pg.y(data), '.-')
ax.set_aspect(2.0)
ert.show(data)

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
# data['rhoa'] = data['r'] * data['k']
ert.show(data)

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
mgr.invert(verbose=True,
           #paraDX=0.3, paraMaxCellSize=10, paraDepth=20, 
           quality=34,
           size=1.,
           lam=20,
           zWeight=0.5, # run inversion: lambda is a smoothness constraint, lower zweight emphasises layering
           robust=False)

#%%
# Now show the result! Use coverage=1 to shut off fading; also, there's probably 
# a matplotlib method I haven't found yet for adding a title. Labelling the colour
# bar does some of that job though.

# plot comparison pseudosections of data & response of inverted model. If inversion
# worked well, they're indistinguishable. Misfit more useful...
mgr.showFit()
data['misfit'] = pg.log(mgr.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, data['misfit'], cMap="bwr", cMin = -100, cMax = 100,
        label='misfit: response/data/err',
        ) #, cMin=-10, cMax=10

# plot inverted section, can change cMap from default 'Spectral_r' (maybe test 
# viridis or similar for accessibility), set colour max/min, probably specify log 
# scale, and toy with coverage?
ax, _ = mgr.showResult(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r", cMin=0.2, cMax=200,
                       label='Apparent Resistivity ohm.m, zwt=0.5, lambda=20, NOT robust')#, coverage=.5) # , cMap="viridis") #, cMap="Spectral_r"); #cMin=1, cMax=500, 
ax.set_ylim(-16,0)
ax.set_xlim(0,140)

# want to set figure size, but ylim and figsize seem to conflict...
# ax.figsize(5,5) # settings for x, y limits and some additional x-axis ticks.

# ax.set_xlim(0,820)
# ax.set_xticks(np.arange(0, 820.1, 50), minor=True) # 50m spaced 'minor' ticks along x

#%%

# A little tester here using a manually created mesh with lines to constrain 
eleclist = list(set(pg.x(data)))
electrodes = np.zeros((len(pg.x(data)),2), dtype=float)
for n,e in enumerate(eleclist):
    electrodes[n] = [e, 0.]#, 0.]
    
spacing = 2.

world = mt.createWorld(start=[min(pg.x(data))-5*spacing, 0.], 
                       end=[max(pg.x(data))+5*spacing, -10*spacing],
                       worldMarker=True)
pg.show(world)

line1 = mt.createLine(start=[min(pg.x(data))-4*spacing, -1],
                                end=[max(pg.x(data))+4*spacing, -1],
                                marker=2) # marker>0 means it functions as a constraint
line2 = mt.createLine(start=[min(pg.x(data))-4*spacing, -13],
                                end=[max(pg.x(data))+4*spacing, -13],
                                marker=3) # marker>0 means it functions as a constraint

world += line1 + line2

pg.show(world)

#%%

for p in electrodes:
    world.createNode(p)
    world.createNode(p -[0, 0.1])

mesh = pg.meshtools.createMesh(world, quality=34, size=1.)

ax, _ = pg.show(mesh)
# ax.set_xlim(-20, 150);
# ax.set_ylim(-100, 10);
#%%
mgrConstrained = ert.ERTManager()
mgrConstrained.invert(data=data,
                      verbose=True, 
                      lam=20, 
                      mesh=mesh, 
                      zWeight=0.5, 
                      robust=False)

#%%
ax, _ = mgrConstrained.showResult(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r", cMin=0.2, cMax=200,
                       label='Apparent Resistivity ohm.m, zwt=0.5, lambda=20, robust');
ax.set_ylim(-16,0)
ax.set_xlim(0,140)


misfit_Const = pg.log(mgrConstrained.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, misfit_Const, cMap="bwr", cMin=-50, cMax=50, label='misfit: response/data/err')

#%%
"""
Now...let's try using the "prior" data...

"""
kw = dict(cMin=0.2, cMax=200, logScale=True, cMap="Spectral_r", 
          xlabel="x (m)", ylabel="z (m)")

x,z,r = prior['x'], prior['z'], prior['r']

ax, cb = mgr.showResult(**kw)
zz = np.abs(z)
iz = np.argsort(z)
dz = np.diff(zz[iz])
thk = np.hstack([dz, dz[-1]])
ztop = -zz[iz[0]] - dz[0]/2
colkw = dict(x=x[0], val=r[iz], thk=thk, width=4, ztopo=ztop)
draw1DColumn(ax, **colkw, **kw)
ax.grid(True)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)

posVec = [pg.Pos(pos) for pos in zip(x, z)]
para = pg.Mesh(mgr.paraDomain)  # make a copy
para.setCellMarkers(pg.IVector(para.cellCount()))
fopDP = PriorModelling(para, posVec)

fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
resSmooth = fopDP(mgr.model)
ax.semilogx(resSmooth, z, label="ERT")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid(True)
ax.legend()

mgr.inv.setRegularization(2, correlationLengths=[40, 4])
mgr.invert()
ax, cb = mgr.showResult(**kw)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)
draw1DColumn(ax, **colkw, **kw)
resGeo = fopDP(mgr.model)

fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
ax.semilogx(resSmooth, z, label="ERT smooth")
#ax.semilogx(res2, z, label="ERT aniso")
ax.semilogx(resGeo, z, label="ERT geostat")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid()
ax.legend()

inv = pg.Inversion(fop=fopDP, verbose=True)
inv.mesh = para
tLog = pg.trans.TransLog()
inv.modelTrans = tLog
inv.dataTrans = tLog
inv.setRegularization(correlationLengths=[40, 4])
model = inv.run(r, relativeError=0.2)
ax, cb = pg.show(para, model, **kw)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)
draw1DColumn(ax, **colkw, **kw)

fopJoint = JointModelling([mgr.fop, fopDP])
# fopJoint.setMesh(para)
fopJoint.setData([data, pg.Vector(r)])  # needs to have .size() attribute! (?)
model = pg.Vector(para.cellCount(), 100)
response = fopJoint(model)
respERT = response[:data.size()]
respDP = response[data.size():]
print(respDP)

fopJoint.createJacobian(model)  # works
J = fopJoint.jacobian()
print(type(J))  # wrong type
ax, cb = pg.show(J)
print(J.mat(0))
ax, cb = pg.show(J.mat(1), markersize=4)

#%%%

dataVec = np.concatenate((data["rhoa"], r))
errorVec = np.concatenate((data["err"], np.ones_like(r)*0.2))
inv = pg.Inversion(fop=fopJoint, verbose=True)
transLog = pg.trans.TransLog()
inv.modelTrans = transLog
inv.dataTrans = transLog
inv.run(dataVec, errorVec, startModel=model)
ax, cb = pg.show(para, inv.model, **kw)
draw1DColumn(ax, **colkw, **kw)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)

inv.setRegularization(2, correlationLengths=[40, 4])
model = inv.run(dataVec, errorVec, startModel=model)
ax, cb = pg.show(para, model, **kw)
draw1DColumn(ax, **colkw, **kw)
ax.set_ylim(-20,0)
ax.set_xlim(0,140)

respERT = inv.response[:data.size()]
misfit = - respERT / data["rhoa"] * 100 + 100
ax, cb = ert.show(data, misfit, cMap="bwr", cMin=-5, cMax=5)

respDP = inv.response[data.size():]
fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
# resMesh = pg.interpolate(srcMesh=para, inVec=inv.model, destPos=posVec)
# ax.semilogx(resMesh, z, label="ERT+DP")
ax.semilogx(respDP, z, label="response")
ax.grid(True)
ax.legend()

"""

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