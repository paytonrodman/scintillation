#!/usr/bin/python

# import all necessary packages
import matplotlib
from matplotlib import pyplot
import json
import mx.DateTime
import numpy as np
import sys

''' Bins input1 into bin defined by the parameters of bin_set, and appends
connected values input2 and input3 into separate lists with matching indices'''
def freq_binning(input1, input2, input3, bin_set, input1_bin, input2_bin, input3_bin):
    #if frequency is less than a bin upper limit, place frequency and
    #corresponding flux, epoch into that bin.
    l = 0
    k = 0
    while l < len(input1):
        if input1[l] < 4.0 or input1[l] > 8.5:
            l += 1
        elif input1[l] < bin_set[k] and input1[l] >= 4.0:
            fr = input1[l]
            input1_bin[k].append(fr)
            fl = input2[l]
            input2_bin[k].append(fl)
            if not input3[-1].isdigit():
                input3 = input3[:-1]
            ep = mx.DateTime.strptime(input3, '%Y-%m-%d').mjd
            input3_bin[k].append(ep)
            l += 1
        else:
            k += 1
    return input1_bin, input2_bin, input3_bin

'''Adds an identifier to the beginning of each bin of input, to retain relationships
between bins'''
def add_ident(input1, ident):
    m = 0
    identifier = ident + repr(m)
    for i in range(len(input1)):
        if ident + repr(m) not in input1[m]:
            input1[m].insert(0, ident + repr(m))
            m += 1
    return input1

'''Within dictionary key 'a', creates a dictionary with key 'type_ident',
then appends 'add_on' under this key'''
def append_dict(sourcename, type_ident, add_on ):
    full_dict[sourcename][type_ident] = []
    full_dict[sourcename][type_ident].append(add_on)





'''Extracts raw data for sample, and outputs a reduced and organized version after
binning'''
infilename = sys.argv[1] #infilename must be manually updated to most recent version
outfilename = sys.argv[2] #extracted and binned values, eventually saved to form: sourcename_[outfilename].json
#read .json data for sources and get handle for each source
data = json.load(open(infilename))
sources = data.keys()
#create empty outfile (json) for saving modulation indices


source_no = 0
for i in sources:
    outfile = []
    outfilenamealt = sources[source_no] + '_' + outfilename + '.json'
    print outfilenamealt
    print 'Source number [', source_no + 1, '] out of [', len(sources), ']'
    ra = data[sources[source_no]]['rightAscension']
    dec = data[sources[source_no]]['declination']
    #get all epochs
    epoch = data[sources[source_no]]['epochs']
    if len(epoch) > 0:
        epoch_mjd = []
        for j in epoch:
            if not j[-1].isdigit():
                j = j[0:-1]
            epoch_mjd.append(mx.DateTime.strptime(j, '%Y-%m-%d').mjd)
        #create set of frequency bins
        start = 4.5
        stop = 8.5
        step = (stop - start)/39
        bin_set = list(np.arange(start, stop + step, step)) #must add 1*(step) to desired end-point.

        freq_list = [[] for _ in range(len(bin_set))]
        flux_list = [[] for _ in range(len(bin_set))]
        epoch_list = [[] for _ in range(len(bin_set))]
        epoch_count_list = [[] for _ in range(len(bin_set))]

        epoch_no = 0
        for i in epoch:
            current_epoch = i
            #get list of all flux and frequency pairs for each epoch
            flux_freq = data[sources[source_no]]['data'][epoch_no]['fluxDensityData'][0]['data']
            flux_freq_transpose = np.transpose(flux_freq)
            #separate pairs of flux and freq into individual lists
            try:
                freq = flux_freq_transpose[0]
                flux = flux_freq_transpose[1]
            except:
                freq = []
                flux = []
            #bin frequencies, fluxes, and epochs
            [freq_list, flux_list, epoch_list] = freq_binning(freq, flux, current_epoch, bin_set, freq_list, flux_list, epoch_list)

            epoch_no += 1
        index = 0
        for i in epoch_list:
            #total number of individual epochs
            epoch_count = len(np.unique(i))
            epoch_count_list[index].append(epoch_count)
            index += 1
        #appends an identifying integer at the beginning of each bin,
        # which later become the keys when list is converted to a dictionary
        freq_list = add_ident(freq_list, 'freq')
        flux_list = add_ident(flux_list, 'flux')
        epoch_list = add_ident(epoch_list, 'epoch')
        epoch_count_list = add_ident(epoch_count_list, 'epochCount')

        #create nested dictionary of data
        #epoch_and_flux = {'epoch': epoch_list,'flux': flux_list}
        pos = {'ra': ra[0], 'dec': dec[0]}
        freq_ident = {k[0]: k[1:] for k in freq_list}
        flux_ident = {k[0]: k[1:] for k in flux_list}
        epoch_ident = {k[0]: k[1:] for k in epoch_list}
        epoch_count_ident = {k[0]: k[1:] for k in epoch_count_list}
        full_dict = {}
        source_name = [sources[source_no]]

        #append all dictionary and list data to the main dictionary
        for a in source_name:
            full_dict[a] = {}
            append_dict(a, 'pos', pos)
            append_dict(a, 'freq', freq_ident)
            append_dict(a, 'flux', flux_ident)
            append_dict(a, 'epoch', epoch_ident)
            append_dict(a, 'epochCount', epoch_count_ident)
            full_dict[a]['maxNumEpochs'] = len(epoch)
            full_dict[a]['epochList'] = epoch_mjd

        #save new data structure
        with open(outfilenamealt, mode='w') as empty:
            entry = full_dict
            outfile.append(entry)
            json.dump(outfile, empty)
   source_no += 1
