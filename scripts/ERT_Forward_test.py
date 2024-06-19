# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:14:38 2024

@author: matthew.noteboom

A fairly simple routine using pyGimli to forward model a layered earth, then
return an inversion result from the forward modelled resistivity.

Might be interesting to look into dipping layers, and polygon targets...

"""

# import useful/necessary packages
import pygimli as pg
import numpy as np
import sys
from pygimli.physics import ert
import pygimli.meshtools as mt

pyversion = sys.version

# Set target folder for input file (if using field data for electrode locations)
# and, if required, for saving forward model output.
datafolder = ".\\exampledata\\Lubango\\"

#%%%
# Load data to container from BERT format to use Acacia Grad electrode config...
# data = ert.load(datafolder+"2024-03-07 OPTHOF-ROLL_GradientXL_2.ohm") 
# elecs=pg.x(data)
# config = 'gr'

# ...or use an automated dip-dip (dd), Wenner (wa), Schlum (slm), P-dip (pd).
# To get a list of configs recognised by createData, input nonsense schemename,
# e.g. 'abc'.

# Currently has electrodes every 5m from 0-100m, in a dipole-dipole scheme
elecs=np.linspace(start=0, stop=100, num=21)
config = 'dd'

scheme = ert.createData(elecs=elecs, schemeName=config, maxSeparation=10)

# Create 'world' with layers, bodies, etc.
# Currently creates a space from -50m to 150m, with depth 100m, with layers 
# defined at 2m and 10m depth above a halfspace
elecrange = np.max(elecs) - np.min(elecs)
spacing = (max(elecs) - min(elecs))/(len(elecs) -1)
world = mt.createWorld(start=[min(elecs)-5*spacing, 0.], 
                       end=[max(elecs)+5*spacing, -10*spacing],
                       layers=[-2., -10.], worldMarker=True)

"""
Examples of bodies from pyGimli notebook:
block = mt.createCircle(pos=[-5, -3.], radius=[4, 1], marker=4,
                        boundaryMarker=10, area=0.1)
poly = mt.createPolygon([(1,-4), (2,-1.5), (4,-2), (5,-2),
                         (8,-3), (5,-3.5), (3,-4.5)], isClosed=True,
                         addNodes=3, interpolate='spline', marker=5)
"""

# Show the 'world' for review
pg.show(world)

#%%
# Create mesh nodes at electrodes, and, from pyGimli recommendations, at 1/10
# the electrode spacing
for p in scheme.sensors():
    world.createNode(p)
    world.createNode(p - [0, 0.1])

# Create mesh for 'world'
mesh = mt.createMesh(world, quality=34)

# Set resistivities for 2 layers and basement halfspace
rhomap = [[1,10.],[2,50.],[3,150.]]

# Take a look at the mesh and the resistivity distribution
ax, _ = pg.show(mesh, data=rhomap, label=pg.unit('resistivity Ohm.m'),
                showMesh=True, cMap="Spectral_r")
ax.set_ylim(-5*spacing,0)

#%%
# Run forward model with inputs of mesh, electrode scheme, resistivities and noise
fwddata = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1,
                    noiseAbs=1e-6, seed=1337)

# set information fields
pg.info(np.linalg.norm(fwddata['err']), np.linalg.norm(fwddata['rhoa']))
pg.info('Simulated data', fwddata)
pg.info('The data contains:', fwddata.dataMap().keys())

pg.info('Simulated rhoa (min/max)', min(fwddata['rhoa']), max(fwddata['rhoa']))
pg.info('Selected data noise %(min/max)', min(fwddata['err'])*100, max(fwddata['err'])*100)

# show forward model result - apparent resistivity pseudosection
# ert.show(data)
ert.show(fwddata)
# ert.show(fwddata, -fwddata['k'], logScale=True)

#%%
# Now invert with ERTManager to compare to input model
# call manager
mgr = ert.ERTManager(fwddata)

# Run inversion (very fast with 'very good' model data and simple geol model)
inv = mgr.invert(lam=20, verbose=True, zWeight=0.1)
# np.testing.assert_approx_equal(mgr.inv.chi2(), 0.7, significant=1)

# plot inversion result with input/output pseudosections
mgr.showResultAndFit()
# plot just inversion result with fixed axis limits
ax, _ = mgr.showResult(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r")
ax.set_xlim(-10,110)
# ax.set_ylim(-25,0)

#%%
# Now try with some lines to constrain inversion...this is quite instructive
# on the influence of even a very simple constraint - a boundary with no fixed
# values either side.

# Create a new 'world', and two layer boundaries (just lines)
geo = mt.createWorld(start=[-50., 0.], end=[150., -100.], worldMarker=True)
line1 = pg.meshtools.createLine(start=[-40, -2], end=[140, -2], marker=2) # marker>0 means it functions as a constraint
line2 = pg.meshtools.createLine(start=[-40, -10], end=[140, -10], marker=3) # marker>0 means it functions as a constraint

# Add the boundaries to the 'world'
geo += line1 + line2

# Show the model without properties
pg.show(geo)

invmesh = pg.meshtools.createMesh(geo, quality=34, size=1.)
ax, _ = pg.show(invmesh)
ax.set_xlim(-50, 150);
ax.set_ylim(-100, 10);

mgrConstrained = ert.ERTManager(fwddata)
coninv = mgrConstrained.invert(verbose=True, lam=20, mesh=invmesh)

mgrConstrained.showResultAndFit()

ax, _ = mgrConstrained.showResult(xlabel="x (m)", ylabel="z (m)", cMap="Spectral_r")
ax.set_xlim(-10,110)
ax.set_ylim(-40,0)
