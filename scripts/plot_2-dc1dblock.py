#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
VES inversion for a blocky model
================================

This tutorial shows how an built-in forward operator is used for inversion.
A DC 1D (VES) modelling is used to generate data, noisify and invert them.
"""

# %%%
# We import numpy, matplotlib and the 1D plotting function
import numpy as np
import pygimli as pg
from pygimli.physics import VESManager

# %%%
## some definitions before (model, data and error)
## simulated ab & mn...
ab2 = np.logspace(0, 3.3, 20)  # AB/2 distance (current electrodes)
mn2 = ab2/3
## ...or examples, in this case from one Puntland line
# ab2 = np.array([1.5, 2.1, 3., 4.2, 6., 9., 13.5, 20., 20., 30., 30., 45., 66.,
#                 100., 150., 150., 220., 220., 330., 500., 750., 1000. ])
# mn2 = np.array([1., 1., 1., 1., 1., 1., 1., 1., 12., 1., 12., 12., 12., 12.,
#                 12., 90., 12., 90., 90., 90., 90., 90 ])

###############################################################################

# %%%
# define a synthetic model and do a forward simulation including noise
# synres = [100., 500., 30., 800.]  # synthetic resistivity
# synthk = [5., 35., 60.]  # synthetic thickness (nlay-th layer is infinite)

synres = [10., 100., 5., 1000.]  # synthetic resistivity
synthk = [20., 50., 100.]  # synthetic thickness (nlay-th layer is infinite)

# synres = [100., 100., 5., 1000.]  # synthetic resistivity
# synthk = [20., 200., 50.]  # synthetic thickness (nlay-th layer is infinite)

###############################################################################
# the forward operator can be called by f.response(model) or simply f(model)
synthModel = synthk + synres  # concatenate thickness and resistivity
ves = VESManager()
rhoa, err = ves.simulate(synthModel, ab2=ab2, mn2=mn2,
                         noiseLevel=0.03, seed=1337)

# %%%
ves.invert(data=rhoa, error=err, ab2=ab2, mn2=mn2,
           nLayers=4,
           # startModel=[ 10., 100., 5., 1000.,
           #              20., 150., 100. ],
           lam=1000, lambdaFactor=0.8
           )

# %%%
# show estimated&synthetic models and data with model response in 2 subplots
fig, ax = pg.plt.subplots(ncols=2, figsize=(8, 6))  # two-column figure
ves.showModel(synthModel, ax=ax[0], label="synth")#, plot="semilogy", zmax=1.5*sum(synthk))
ves.showModel(ves.model, ax=ax[0], label="model")#, zmax=1.5*sum(synthk))
ves.showData(rhoa, ax=ax[1], label="data", color="C0", marker="x")
ves.showData(ves.inv.response, ax=ax[1], label="response", color="C1")

ax[0].set(ylim=(1000,1)) # For model subplot; x is res, y is depth (xlim=(1,10000), )
ax[1].set(ylim=(2000,1))  # For data/response subplot, x is res, y is ab/2 (xlim=(1,1000), )
