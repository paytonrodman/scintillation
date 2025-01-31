#!/usr/bin/python

# import all necessary packages
import matplotlib
from matplotlib import pyplot
import scipy.signal
import json
import numpy as np
import sys
import statsmodels.robust

infilename = sys.argv[1] #infilename should be the outfile from 'calculate_sf.py'. The contents of this file are altered, to remove outliers
outfilename = sys.argv[2] #outfile lists any source which had outliers removed, and the frequency bins of the outliers
outfilename = outfilename + '.json'
data = []
outfile = []
full = {}

with open(infilename) as f:
    for line in f:
        data.append(json.loads(line))

for source_no in range(len(data)):
    sourcename = data[source_no][0].keys()
    sourcename = sourcename[0]
    print 'Source number [', source_no + 1, '] out of [', len(data), ']'

    data_type = data[source_no][0][sourcename] #gives dictionary to 'tau', 'sf', or 'mod'
    mod = data_type['mod']
    tau = data_type['tau']
    sf = data_type['sf']
    mod_numel = data_type['mod_numel']
    num_epochs = data_type['num_epochs']
    pos = data_type['pos']
    mod_list = []

    for i in range(len(mod[0])):
        string = 'mod' + repr(i)
        mod_list.append(mod[0][string])

    #replace NaN values in mod_list with its nearest lefthand neighbour
    lastgood = np.nanmean(mod_list)

    for i,value in enumerate(mod_list):
        if not np.isnan(value):
            lastgood = value
        else:
            mod_list[i] = lastgood
            
    #apply median window filter and median absolute deviation
    filt = scipy.signal.medfilt(mod_list, 15)
    sigma = statsmodels.robust.mad(mod_list)
    flagged = (np.abs(mod_list - filt) / sigma) > 5

    #if there are any outliers flagged
    if True in flagged:
        flag_idx = [i for i,x in enumerate(flagged) if x == True]
        for j in flag_idx:
            mod_string = 'mod' + repr(j)
            sf_string = 'sf' + repr(j)
            tau_string = 'tau' + repr(j)
            numel_string = 'numel' + repr(j)
            data[source_no][0][sourcename]['sf'][0][j][sf_string] = []
            data[source_no][0][sourcename]['tau'][0][j][tau_string] = []
            data[source_no][0][sourcename]['mod'][0][mod_string] = np.NaN
            data[source_no][0][sourcename]['mod_numel'][0][numel_string] = 0
        full[sourcename] = flag_idx

with open(outfilename, mode='w') as f:
    outfile.append(full)
    json.dump(outfile, f)
with open(infilename, mode='w') as f:
    f.write(json.dumps(data))
