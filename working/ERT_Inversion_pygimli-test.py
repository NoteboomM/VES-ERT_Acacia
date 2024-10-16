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


"""
# import useful/necessary packages
import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from pygimli.physics import ert
# import pybert as pb

pyversion = sys.version

# This version works with a temporary testing folder, but can be changed to any path
# (although eventually there is probably a length limit).
os.chdir("C:\\temp\\test-ERT2D")

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
data['rhoa'] = data['r'] * data['k']
ert.show(data)

#%% OPTIONAL BLOCK!
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

#%%

# data['err'] = ert.estimateError(data, # simulate an error field if necessary
#                                 absoluteUError=0.00005, # 50µV
#                                 relativeError=0.03)  # 3% 

# plot a pseudosection of the error
ert.show(data, data['err']*100, label="error [%]")
data.remove(data["err"] >= 5/100 ) # filter extremes?

#%%
# call the ERTManager class
mgr = ert.ERTManager(data)
# run the inversion, with some settings for mesh, lambda, smoothness
mgr.invert(verbose=True,
           #paraDX=0.3, paraMaxCellSize=10, paraDepth=20, 
           quality=34,
           size=1.,
           lam=20,
           zWeight=1, # run inversion: lambda is a smoothness constraint, lower zweight emphasises layering
           robust=False)

#%%
# Now show the result! Use coverage=1 to shut off fading; also, there's probably 
# a matplotlib method I haven't found yet for adding a title. Labelling the colour
# bar does some of that job though.

# plot comparison pseudosections of data & response of inverted model. If inversion
# worked well, they're indistinguishable. Misfit more useful...
mgr.showFit()
data['misfit'] = pg.log(mgr.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, data['misfit'], cMap="bwr", label='misfit: response/data/err') #, cMin=-10, cMax=10

# plot inverted section, can change cMap from default 'Spectral_r' (maybe test 
# viridis or similar for accessibility), set colour max/min, probably specify log 
# scale, and toy with coverage?
ax, _ = mgr.showResult(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r", cMin=0.4, cMax=1000,
                       label='Apparent Resistivity ohm.m, lambda=20, NOT robust')#, coverage=.5) # , cMap="viridis") #, cMap="Spectral_r"); #cMin=1, cMax=500, 
ax.set_ylim(-25,0)
ax.set_xlim(0,140)

# want to set figure size, but ylim and figsize seem to conflict...
# ax.figsize(5,5) # settings for x, y limits and some additional x-axis ticks.

# ax.set_xlim(0,820)
# ax.set_xticks(np.arange(0, 820.1, 50), minor=True) # 50m spaced 'minor' ticks along x

#%%

# A little tester here using a manually created mesh with a line to constrain between
# higher resistivity ~0-10m depth and lower resistivity at greater depth with
# 2024-03-07 test data from Acacia Terrameter
geo = pg.meshtools.createParaMeshPLC(data, paraMaxCellSize=100)
line = pg.meshtools.createLine(start=[30, -10], end=[145, -10], marker=1) # marker>0 means it functions as a constraint

ax, _ = pg.show(geo);
ax.set_xlim(-20, 150);
ax.set_ylim(-50, 10);

geo += line
mesh = pg.meshtools.createMesh(geo, quality=34, size=1.)
ax, _ = pg.show(mesh)
ax.set_xlim(-20, 150);
ax.set_ylim(-100, 10);

mgrConstrained = ert.ERTManager()
mgrConstrained.invert(data=data, verbose=True, lam=200, mesh=mesh, zWeight=1, robust=True)

ax, _ = mgrConstrained.showResult(cMap="Spectral_r", cMin=2, cMax=200, xlabel="x (m)", ylabel="z (m)",
                                  label='Apparent Resistivity ohm.m, lambda=200, robust');
ax.set_ylim(-25,0)
ax.set_xlim(0,140)


misfit_Const = pg.log(mgrConstrained.inv.response / mgr.data["rhoa"]) / data["err"]
pg.show(data, misfit_Const, cMap="bwr", cMin=-50, cMax=50, label='misfit: response/data/err')
