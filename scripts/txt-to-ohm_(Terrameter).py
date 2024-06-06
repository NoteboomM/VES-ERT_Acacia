# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 09:53:18 2024

@author: matthew.noteboom

Routine to convert stg files to .ohm files (from the BERT working group:
http://resistivity.net/bert/data_format.html ) for the Lubango survey. Should be 
able to operate with minimal changes for any stg data files.

"""
#%%
import os,sys
# store version for reference (habit - sometimes useful for troubleshooting)
pyversion = sys.version

# change to working directory
os.chdir("C:/temp/test-ERT2D")

#%%
# open and load stg file; if required later can maybe improve with a prompt or listing?
# with open("DD4384(1-6).stg", 'r') as stgfile:
with open("2024-03-07 OPTHOF-ROLL_GradientXL_2.txt", 'r') as lsfile:
    lsdata = lsfile.readlines()[1:]
# open and load terrain file. Left this in from another script in case we have 
# something on Terrameter surveys, but it's probably not very often important in NL!
# terrain can influence geometric factor(s) relative to simple analytic calc.
# with open("ALT#5-1.trn", 'r') as trnfile:
#     trndata = trnfile.readlines()[3:]

# load terrain data into dictionary to later call to build electrode location block of output
# trndict = {}
# for line in trndata:
#     x = line.split()[0]
#     z = line.split()[1]
#     trndict[x] = z

# Simple approach for the number of measurements 
numstations = len(lsdata)

# initialise output data block with number of stations and channel labels.
# Leaving old structure in place to handle IP if/when it happens
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
    lineout=(f'{eleclist.index(a):<5}{eleclist.index(b):<5}'
             f'{eleclist.index(m):<5}{eleclist.index(n):<5}'
             f'{r:<14.5E}{err:<14.5E}{rhoa:<.5E}\n')
    
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
with open("2024-03-07 OPTHOF-ROLL_GradientXL_2.ohm", 'w') as outfile:
    outfile.write(outdata)
