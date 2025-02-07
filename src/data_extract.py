#!/usr/bin/python
#
# data_extract.py
#
# A program to read raw ATESE data and create data files containing the relevant
# variables.
#
# Usage: python data_extract.py -i [input_file] -o [output_dir]
#
# Output:
#    - a .json file for each source including the source name, location, frequency
#      and flux data per epoch, and details about each epoch
#

# Python standard modules
import argparse
import sys

# Other Python modules
import json
import mx.DateTime as DateTime
import numpy as np
from funcs import add_ident, append_dict, bin_frequency

def main(**kwargs):
    # Read .json data for sources and get handle for each source
    data = json.load(open(kwargs['input_file']))
    sources = data.keys()

    source_no = 0
    for i in sources:
        outdata = []

        # Get source information (RA, Dec, epochs)
        ra = data[sources[source_no]]['rightAscension']
        dec = data[sources[source_no]]['declination']
        epoch = data[sources[source_no]]['epochs']

        if len(epoch) > 0:
            # Convert epochs to mjd date format
            epoch_mjd = []
            for j in epoch:
                if not j[-1].isdigit():
                    j = j[0:-1]
                epoch_mjd.append(DateTime.strptime(j, '%Y-%m-%d').mjd)

            # Create set of frequency bins
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

                # Get list of all flux and frequency pairs for each epoch
                flux_freq = data[sources[source_no]]['data'][epoch_no]['fluxDensityData'][0]['data']
                flux_freq_transpose = np.transpose(flux_freq)

                # Separate pairs into individual lists
                try:
                    freq = flux_freq_transpose[0]
                    flux = flux_freq_transpose[1]
                except:
                    freq = []
                    flux = []

                # Bin frequencies, fluxes, and epochs
                [freq_list, flux_list, epoch_list] = bin_frequency(freq, flux, current_epoch, bin_set, freq_list, flux_list, epoch_list)

                epoch_no += 1

            # Get number of individual epochs
            index = 0
            for i in epoch_list:
                epoch_count = len(np.unique(i))
                epoch_count_list[index].append(epoch_count)
                index += 1

            # Append an identifying integer at the beginning of each bin,
            # which later become the keys when list is converted to a dictionary
            freq_list = add_ident(freq_list, 'freq')
            flux_list = add_ident(flux_list, 'flux')
            epoch_list = add_ident(epoch_list, 'epoch')
            epoch_count_list = add_ident(epoch_count_list, 'epochCount')

            # Create nested dictionary of data
            # epoch_and_flux = {'epoch': epoch_list,'flux': flux_list}
            pos = {'ra': ra[0], 'dec': dec[0]}
            freq_ident = {k[0]: k[1:] for k in freq_list}
            flux_ident = {k[0]: k[1:] for k in flux_list}
            epoch_ident = {k[0]: k[1:] for k in epoch_list}
            epoch_count_ident = {k[0]: k[1:] for k in epoch_count_list}
            full_dict = {}
            source_name = [sources[source_no]]

            # Append all dictionary and list data to the main dictionary
            for sname in source_name:
                full_dict[sname] = {}
                full_dict = append_dict(full_dict, sname, 'pos', pos)
                full_dict = append_dict(full_dict, sname, 'freq', freq_ident)
                full_dict = append_dict(full_dict, sname, 'flux', flux_ident)
                full_dict = append_dict(full_dict, sname, 'epoch', epoch_ident)
                full_dict = append_dict(full_dict, sname, 'epochCount', epoch_count_ident)
                full_dict[sname]['maxNumEpochs'] = len(epoch)
                full_dict[sname]['epochList'] = epoch_mjd


            # Save new data structure
            savename = kwargs['output_dir'] + '/' + sources[source_no] + '.json'
            with open(savename, mode='w') as file:
                outdata.append(full_dict)
                json.dump(outdata, file)

       source_no += 1


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract raw data for sample,\
                      and output a reduced and organized version after binning')
    parser.add_argument('-i', '--input_file',
                        type=argparse.FileType('r'),
                        nargs=1,
                        default=None,
                        help='Raw data file from ATESE pipeline')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        default=None,
                        help='Extracted and binned values')
    args = parser.parse_args()

    main(**vars(args))
