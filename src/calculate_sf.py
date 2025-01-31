#!/usr/bin/python

# import all necessary packages
import matplotlib
from matplotlib import pyplot
import json
import mx.DateTime
import numpy as np
import sys

''' Creates a string for identifying the bin number (index) of name1'''
def bin_str(name1, index):
    name1_str = name1 + repr(index)
    return name1_str

''' Bins input1 into bin defined by the parameters of bin_set, and appends
connected values of input2 separate list with matching indices'''
def binning(input1, input2, bin_set, input1_bin, input2_bin):
    #if tau is less than first bin limit, place in first bin otherwise, test if
    #tau is between two bins, and place in corresponding bin
    l = 0
    for i in input1:
        if input1[l] <= bin_set[0]:
            dt = input1[l]
            input1_bin[0].append(dt)
            sf = input2[l]
            input2_bin[0].append(sf)
        else:
            for k in range(1, len(bin_set)):
                if input1[l] <= bin_set[k] and input1[l] > bin_set[k - 1]:
                    dt = input1[l]
                    input1_bin[k].append(dt)
                    sf = input2[l]
                    input2_bin[k].append(sf)
        l += 1
    return input1_bin, input2_bin





'''Calculates the full structure function as a function of time difference (tau),
along with modulation index calculated as (RMS flux / mean flux), which are then
appended to the ouput file'''
infilename = sys.argv[1] #infilename should be the output of 'data_extract.py' - individual files
print 'file in is:', infilename
outfilename = sys.argv[2] #sf/tau calculated values, saved to form: [outfilename].json
# read .json data for sources and get handle for each source
data = json.load(open(infilename))
part = infilename.split("/")[3]
number = part.split("_")[0]
outfilenamealt = number + '_' + outfilename + '.json'

#create empty dictionaries for saving data later
full = {}
#loop through sources
for source_no in range(len(data)):
    outfile = []
    #print 'source number: [', (source_no+1), '] out of ', len(data)
    #get source name and assign a sub-dictionary to it
    sources = data[source_no].keys()
    full[sources[0]] = {}
    sourcename = sources[0]
    if sourcename in sources: #option to select particular sources to calculate
        print 'source is:', sourcename
        #reads data for one freq bin and the assocated flux, epoch bin
        data2 = data[source_no][sourcename]['freq'][0]
        subtau = []
        subsf = []
        submod = {}
        submodel = {}
        subvar = {}
        max_tau_list = []
        for i in range(0, len(data2)):
            epoch_str = bin_str('epoch', i)
            epoch_bin = data[source_no][sourcename]['epoch'][0][epoch_str]
            try:
                max_ep = max(epoch_bin)
                min_ep = min(epoch_bin)
                max_tau = max_ep - min_ep
            except:
                'fail try1'
                max_tau = 0
            max_tau_list.append(max_tau)
        max_tau_all = max(max_tau_list)
        print 'max tau is:', max_tau_all

        #loop through frequency bins
        for ident in range(0, len(data2)):
            #print 'frequency number is: [', (ident+1), '] out of ', len(data2)
            #creates key name for each frequency bin, for calling
            freq_str = bin_str('freq', ident)
            flux_str = bin_str('flux', ident)
            epoch_str = bin_str('epoch', ident)
            mod_str = bin_str('mod', ident)
            numel_str = bin_str('numel', ident)
            var_str = bin_str('var', ident)
            #calls all data for one bin number
            freq_bin = data[source_no][sourcename]['freq'][0][freq_str]
            flux_bin = data[source_no][sourcename]['flux'][0][flux_str]
            epoch_bin = data[source_no][sourcename]['epoch'][0][epoch_str]
            num_epochs = data[source_no][sourcename]['maxNumEpochs']
            pos = data[source_no][sourcename]['pos']
            #calculates the unbinned structure function and corresponding tau values
            #assigns the same epoch to each of the previous slices made
            tau = []
            sf_list = []
            for j in range(0, len(flux_bin)):
                for k in range(j + 1, len(flux_bin)):
                    flux_diff = abs(flux_bin[j] - flux_bin[k])
                    sf = flux_diff**2
                    sf_list.append(sf)
                    time_diff = abs(epoch_bin[j] - epoch_bin[k])
                    tau.append(time_diff)
            try:
                #create bin set for time difference, tau
                start = 0
                stop = max_tau
                step = abs((start - stop)/19)
                #add (step) onto ending point to ensure all points are counted
                try:
                    bin_set = list(np.arange(start, stop + step, step))
                except ZeroDivisionError as error:
                    print 'Only 1 epoch'
                    bin_set = []
                #empty lists for storing the binned values
                tau_binned = [[] for _ in range(len(bin_set))]
                sf_binned = [[] for _ in range(len(bin_set))]
                #bin tau and the sf
                #print len(tau)
                #print len(sf_list)
                [tau_binned, sf_binned] = binning(tau, sf_list, bin_set, tau_binned, sf_binned)
                tau_val_bin = []
                sf_mean_bin = []
                #average tau and sf over each frequency bin, and append to dictionary
                for i in bin_set:
                    halfway  = i + (1/2)*step
                    tau_val_bin.append(halfway)
                tau_str = bin_str('tau', ident)
                subtau_add = {tau_str: tau_val_bin}
                subtau.append(subtau_add)
                for i in sf_binned:
                    if len(i) > 0:
                        av_sf = np.mean(i)
                    else:
                        av_sf = np.NaN
                    sf_mean_bin.append(av_sf)
                sf_str = bin_str('sf', ident)
                subsf_add = {sf_str: sf_mean_bin}
                subsf.append(subsf_add)
            except RuntimeWarning:
                tau_str = bin_str('tau', ident)
                subtau_add = {tau_str: []}
                subtau.append(subtau_add)
                sf_str = bin_str('sf', ident)
                subsf_add = {sf_str: []}
                subsf.append(subsf_add)
            #try calculating the modulation index, with NaN for empty frequency bins
            #will still return warnings
            if len(epoch_bin) > 0:
                numel = len(np.unique(epoch_bin))
                rms_var = np.std(flux_bin, ddof=1)
                mean = np.mean(flux_bin)
                mod_index = rms_var/mean
            else:
                numel = 0
                mod_index = np.NaN
                print 'cannot calculate mod idx for freq bin ', ident
            #create dictionary structure of modulation index data, linked to
            # corresponding frequency bin
            submod_add = {mod_str: mod_index}
            submod.update(submod_add)
            submodel_add = {numel_str: numel}
            submodel.update(submodel_add)
            subvar_add = {var_str: rms_var}
            subvar.update(subvar_add)
        #append dictionary of modulation index data to each source
        for i in sources:
            full[i]['sf'] = []
            full[i]['sf'].append(subsf)
            full[i]['tau'] = []
            full[i]['tau'].append(subtau)
            full[i]['mod'] = []
            full[i]['mod'].append(submod)
            full[i]['mod_numel'] = []
            full[i]['mod_numel'].append(submodel)
            full[i]['rms_var'] = []
            full[i]['rms_var'].append(subvar)
            full[i]['num_epochs'] = num_epochs
            full[i]['pos'] = pos
    #appends sf and modulation index data to the empty file
    with open(outfilenamealt, mode='w') as empty:
        entry = full
        outfile.append(entry)
        json.dump(outfile, empty)
