#!/usr/bin/python

# import all necessary packages
import matplotlib
from matplotlib import pyplot
import json
import mx.DateTime
import numpy as np
import sys

infilename1 = sys.argv[1] #should be the outfile from 'calculate_sf.py' -- all files
infilename2 = sys.argv[2] #should be the directory to the relevant outfiles from data_extract.sf -- individual files

sfdata = []
with open(infilename1) as f:
    for line in f:
        sfdata.append(json.loads(line))

for source_no in range(len(sfdata[0])):
    subrate = {}
    subsfnorm = {}
    sourcename = sfdata[0][source_no][0].keys()
    sourcename = sourcename[0]
    print 'Source number [', source_no + 1, '] out of [', len(sfdata[0]), '] in sfdata'
    sf = sfdata[0][source_no][0][sourcename]['sf']

    fluxdata = json.load(open(infilename2 + sourcename + '_data.json'))
    flux_dict = fluxdata[0][sourcename]['flux'][0]
    flux_all = []
    for freq_bin in flux_dict:
        flux_all.append(flux_dict[freq_bin])
    av_sourceflux = np.mean(flux_all[0])
    sfdata[0][source_no][0][sourcename]['av_flux'] = av_sourceflux
    rms = sfdata[0][source_no][0][sourcename]['rms_var']

    for freq_num in range(len(sf[0])):
        sf_key = 'sf' + repr(freq_num)
        sf_norm_key = 'sfnorm' + repr(freq_num)
        rate_key = 'rate' + repr(freq_num)
        rms_key = 'var' + repr(freq_num)
        sf_bin = sf[0][freq_num][sf_key]
        rms_bin = rms[0][rms_key]

        try:
            rate_bin_norm = [(0.02*(np.sqrt(x/(rms_bin**2)))) for x in sf_bin]
        except ZeroDivisionError as error:
            rate_bin_norm = [np.NaN for x in sf_bin]
        sf_bin_norm = [(np.sqrt(x)/av_sourceflux) for x in sf_bin]

        subrate_add = {rate_key: rate_bin_norm}
        subrate.update(subrate_add)
        subsfnorm_add = {sf_norm_key: sf_bin_norm}
        subsfnorm.update(subsfnorm_add)
    sfdata[0][source_no][0][sourcename]['rate'] = []
    sfdata[0][source_no][0][sourcename]['rate'].append(subrate)
    sfdata[0][source_no][0][sourcename]['sf_norm'] = []
    sfdata[0][source_no][0][sourcename]['sf_norm'].append(subsfnorm)

with open(infilename1, mode='w') as f:
    f.write(json.dumps(sfdata))
