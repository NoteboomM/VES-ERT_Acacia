# -*- coding: utf-8 -*-
"""

VES inversion for a blocky model
================================

Code to invert VES data from Cunene, Angola. Extract data from *.stg files, then
process based on tutorial script from pyGimli that demonstrated fwd & inverse modelling.

Generates short output csv (if required) with resistivity/thickness/depth error,
and 2 subplots with inversion result left, and observation & model response right.

***              Currently set up with my full path                ***
*** edit working directory path to be relative or replace with own ***

Matt Noteboom Dec/2023

"""

# %%% import packages
# We import numpy, pygimli for processing, and os,sys for admin/navigation
import numpy as np
import pygimli as pg
from pygimli.physics import VESManager
import os, sys
from glob import glob
pyversion = sys.version

# %%% vesinvert function
"""
Define vesinvert function - takes arrays of mn/2, ab/2, res, error,
number of layers, a starting model (can be empty), and 'sitecode' for name of
location used for filename & title.
"""

def vesinvert(mn2, ab2, rho, err, layers, startmod, sitecode):
    # build filename stem
    outname = sitecode.split(',')[0] + '_' + str(layers) + 'layers_pygimli-inversion'

    ves = VESManager()
    # run inversion, using startmod if it's been populated.
    if len(startmod) > 0:
        ves.invert(data=rho, err=err, ab2=ab2, mn2=mn2,
                   nLayers=layers,
                   startModel=startmod,
                   lam=1000, lambdaFactor=0.8
                   )
        outname = outname + "-SEGEO"
    else:
        ves.invert(data=rho, err=err, ab2=ab2, mn2=mn2,
                   nLayers=layers,
                   # startModel=startmod,
                   lam=1000, lambdaFactor=0.8
                   )


    # Initialise output string and build lists for output
    outdata = ""
    invthk  = ["Thickness"] + list(ves.model[:(layers-1)]) + [-99.]
    invres  = ["Resistivity"] + list(ves.model[layers-1:])
    deplist = ["Depth"]
    depth   = 0.
    
    # maxres was created for a check of the highest resistivity returned on a 
    # series of inversions. Not essential!
    global maxres
    maxres = max(list(ves.model[layers-1:]))
    
    # calculate RMS and RMS% from observations and forward-modelled response of layered-earth model
    # rmse = np.sqrt(np.mean(np.square(ves.inv.response - rho)))
    rmspe = 100 * np.sqrt(np.mean(np.square(ves.inv.response - rho)))/np.mean(ves.inv.response)

    # print(rmse)
    # create list of depth values for output by summing thicknesses
    for i in range(1,len(ves.model[layers-1:])):
        depth = depth + invthk[i]
        deplist.append(depth)
    deplist = deplist + [-99.]
    
    # Iterate through lists to build up comma-delimited output
    for n in range(0,layers+1):
        outdata = outdata + str(invres[n]) + ',' + str(invthk[n]) + ',' + str(deplist[n]) + '\n'
    outdata = outdata + str(layers) + " Layers, RMS% =, " + str(rmspe)
    # print(outdata)
    # for n in range(len(ab2list)):
    #     print(VESab2[n],ves.inv.response[n])

    # Initialise plots; model on left (ax[0]), observation & response on the right (ax[1])
    fig, ax = pg.plt.subplots(ncols=2, figsize=(8, 6))  # two-column figure, shared y-axis (, sharey=True)
    fig.suptitle(sitecode) # Set title
    
    # Plot model & data
    ves.showModel(ves.model, ax=ax[0], label="model", color="C0", plot="loglog")#, plot="semilogy") # , zmax=1000
    ves.showData(rho, ax=ax[1], label="data", color="C0", marker="x")
    ves.showData(ves.inv.response, ax=ax[1], label="response", color="C1")
    
    # set axis limits to produce plots with consistent scale across sites/areas
    ax[0].set(xlim=(1,10000), ylim=(1250,1)) # x is resistivity, y is depth
    ax[1].set(xlim=(1,1200), ylim=(1250,1))  # x is resistivity, y is ab/2
    
    # save plots as PNG, model summary as CSV. Use comments or triple quite below
    # to switch on/off. Plots will always display inline or in plot panel (Spyder)
    fig.savefig(outname+'.png', dpi=200)
    with open(outname+'.csv', 'w') as outfile:
        outfile.write(outdata)
"""
"""

# %%% load data and call invert function
"""
Set working directory, load files, extract inputs for inversion and convert to arrays.

"""
# set data/output folder
os.chdir("C:/temp")
# look for observation files
filelist = glob("V*.stg")

# part of a test for maximum resistivity across series of locations; not essential
allmax = 0.

# Create empty startmod, which will trigger automatic inversion, otherwise...
startmod = []
# ...replace empty startmodel with list of resistivities and thicknesses as a
# starting model. List should have n resistivities and n-1 thicknesses (bottom
# layer is a half space).
# )
startmod = [ 50., 250., 50., 150., 400., 100., 1250.,
             7.5,  15., 30.,  50.,  90., 250.] # VES02
startmod = [ 4.,   1.,  2., 220., 1100.,
             20,  12., 30.,  80.] # VES03 - bad!
startmod = [ 15.7,   23., 184.1, 1980.1,
            29.38, 96.52, 122.9]  #  VES04
startmod = [ 57.4,    37.,   49.4,   11.9,  136.5, 1813.7,
              7.1,    20.,  32.99, 112.34, 156.12, ] # VES05
startmod = [40.1, 18.7,   8.9,  5.1,  24.5, 112., 1478.3,
             3.1,  2.1, 31.19, 92.8, 60.81, 171.55] #VES06
startmod = [69.3, 41.9, 3.7, 8.7, 2.1, 6.5, 185.2, 3589.2,
            4., 3.5, 11., 22., 56.15, 68., 179. ]  #  VES07

# define number of layers for inversion function based on startmod above
if len(startmod) > 0:
    # if starting model defined:
    layers = int((len(startmod) + 1)/2)
else:
    # else set number of layers required/desired
    layers = 7

# work through list of stg files, or select sample by subscripting list
# load data to arrays, then call inversion function above
for stgfile in filelist[5:6]:#
    print("\n***********************\n", stgfile, "\n***********************\n")
    # create input lists
    rholist = []
    ab2list = []
    mn2list = []
    
    # load file
    with open(stgfile,"r") as infile:
        indata = infile.readlines()

    # break up input data and extract geometry and resistivity values
    for line in indata[3:]:
        lineaslist = [x.strip() for x in line.split(',')]
        # n = lineaslist[0]
        date = lineaslist[2]
        rho = float(lineaslist[7])
        if rho < 0: rho = rho * (-1) # for negative values in VES08 file
        rholist.append(rho)
        site = lineaslist[8]
        ab2neg = float(lineaslist[9])
        ab2 = float(lineaslist[12])
        ab2list.append(ab2)
        mn2neg = float(lineaslist[15])
        mn2 = float(lineaslist[18])
        mn2list.append(mn2)
        # integrity check for geometry figures in file; don't think it's required really
        if (ab2 + ab2neg + mn2 + mn2neg) != 0: 
            print("\n*** Problem with AB and MN values in", stgfile,"\n")
            sys.exit()
    
    if len(startmod) > 0:
        sitecode = site + ", " + date[:4] + "-" + date[4:6] + "-" + date[6:8] + ", SEGEO start model"
    else:
        sitecode = site + ", " + date[:4] + "-" + date[4:6] + "-" + date[6:8] + ", auto start model"

    # convert lists to arrays for VESManager
    VESrho = np.array(rholist)
    VESab2 = np.array(ab2list)
    VESmn2 = np.array(mn2list)
    err = np.ones(len(rholist)) * 0.03 # make obs error array; 3% carried over from pygimli tutorial

    # call inversion function
    vesinvert(VESmn2, VESab2, VESrho, err, layers, startmod, sitecode)
    if maxres > allmax: allmax = maxres
    
