# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 09:53:18 2024

@author: matthew.noteboom
Last update: 2024/09/26

Routine to convert Terrameter txt file(s) to .ohm files (from the BERT working group:
http://resistivity.net/bert/data_format.html ) for processing with other
pyGIMLi scripts & notebooks. Limited testing, but as long as Terrameter outputs
are fairly consistent format from one project to another, it should be stable.

"""
#%%
import os,sys
## store version for reference (habit - sometimes useful for troubleshooting)
pyversion = sys.version

## change to working directory
os.chdir("C:/Temp/211269_OptHof_MSNTesting/2024-03-07 OPTHOF-ROLL_LS214100268_19-21-11")

#%%
# open and load txt file; if required later can maybe improve with a prompt or listing?

with open("2024-03-07 OPTHOF-ROLL_GradientXL_2_edited.txt", 'r', errors='ignore') as lsfile:
    # Check the line index number in brackets in next line. This file has data from line 2,
    # so we load 1 onwards (first line is zero). But other example file from Op 't Hof 2024
    # only has data from line 224 onwards, so we would need lsdata = lsfile.readlines()[223:]
    lsdata = lsfile.readlines()[1:]

## open and load terrain file. I left this in from another script in case we have 
## something on Terrameter surveys, but it's probably not very often important in NL!
## Terrain can influence geometric factor(s) relative to simple analytic calculation
## but the effect is subtle on flat/gentle terrain

# with open("ALT#5-1.trn", 'r') as trnfile:
#     trndata = trnfile.readlines()[3:]

## load terrain data into dictionary to later call to build electrode location block of output
# trndict = {}
# for line in trndata:
#     x = line.split()[0]
#     z = line.split()[1]
#     trndict[x] = z

## Simple approach for the number of measurements, assuming no 'footer' in the file
## after the recorded data.
numstations = len(lsdata)

## initialise output data block with number of stations and channel labels.
## Leaving old structure in place to handle IP if/when it happens
# if "IP" in stgdata[0]:
#     datablock = str(numstations) + "\n# a b m n r ip/ms err/%\n"
# then for files without IP data:
# else:
datablock = str(numstations) + "\n# a b m n r err/% rhoa\n"
    
# create an empty list for the electrode locations
eleclist = []

# loop through lines in input data purely to extract and list electrode locations
for line in lsdata:
    # if 'USER' in line: # data lines only
    # load in electrode locations
    a = float(line.split('\t')[4])
    b = float(line.split('\t')[7])
    m = float(line.split('\t')[10])
    n = float(line.split('\t')[13])
    
    # populate list of electrode locations with multiple 'if's as we only want to list once
    if a not in eleclist:
        eleclist.append(a)
    if b not in eleclist:
        eleclist.append(b)
    if m not in eleclist:
        eleclist.append(m)
    if n not in eleclist:
        eleclist.append(n)

# sort the list of electrode locations and count (for header, and QC)
eleclist.sort()
numelecs = len(eleclist)

# loop through input data again for formatted data
for line in lsdata:
    # if 'USER' in line: # data lines only
    
    # load in electrode locations
    a = float(line.split('\t')[4])
    b = float(line.split('\t')[7])
    m = float(line.split('\t')[10])
    n = float(line.split('\t')[13])
    
    # load in most physical data and error
    r = float(line.split('\t')[24])
    i = float(line.split('\t')[20]) # in mA
    u = float(line.split('\t')[22]) # in V
    rhoa = float(line.split('\t')[26])
    err = float(line.split('\t')[25])
    
    # load IP if present, and format outputs
    # if "IP" in line:
    #     ip = float(line.split(',')[30])
    
    #     # format parameters for data block in output file with IP
    #     lineout=(f'{eleclist.index(a):<5}{eleclist.index(b):<5}'
    #              f'{eleclist.index(m):<5}{eleclist.index(n):<5}'
    #              f'{r:<14.5E}{ip:<14.5E}{err:<.5E}\n')
    # else:
    # format parameters for data block in output file without IP

    lineout=(f'{eleclist.index(a)+1:<5}{eleclist.index(b)+1:<5}'
             f'{eleclist.index(m)+1:<5}{eleclist.index(n)+1:<5}'
             f'{r:<14.5E}{err:<14.5E}{rhoa:<.5E}\n')
    # Note to self: need to add 1 to each index so that electrodes are not 
    # zero-indexed, as 0 electrode position is treated as infinite
    
    # append to data block
    datablock += lineout
        
# build electrode location block for output: first header
elecblock = str(numelecs) + "\n# x z\n"
# then list of electrodes by x location, with z information from terrain dictionary
for x in eleclist:
    # elecblock += f'{int(x):<5}{trndict[str(x)]:<5}\n'
    elecblock += f'{int(x):<5}{"0.0":<5}\n'

# concatenate electrode list block and data block
outdata = elecblock + datablock

# write to output file
with open("2024-03-07 OPTHOF-ROLL_GradientXL_2_edited.ohm", 'w') as outfile:
    outfile.write(outdata)
