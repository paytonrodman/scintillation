#!/usr/bin/python
#
# calculate_sf.py
#
# A program to calculate the Structure Function as a function of time difference
# (tau), along with the modulation index.
#
# Usage: python calculate_sf.py -i [input_dir] -o [output_dir]
#
# Output:
#    - a .json file with the Structure Function and modulation index of each source
#

# Python standard modules
import argparse
import sys
import glob
import warnings

# Other Python modules
import json
import numpy as np
from funcs import binning2, make_index_str


def main(**kwargs):
    json_files = glob.glob(kwargs['input_dir'] + '/*.json')
    for jfile in json_files:
        # Read .json data for sources and get handle for each source
        jdata = json.load(open(jfile))
        part = jfile.split("/")[3]
        number = part.split("_")[0]
        outname = number + '_SF.json'

        # Create empty dictionaries for saving data later
        full = {}

        # Loop through sources
        for source_no in range(len(jdata)):
            outdata = []

            # Get source name and assign a sub-dictionary to it
            source = jdata[source_no].keys()
            full[source[0]] = {}

            # Read data for one freq bin and the assocated (flux, epoch) bin
            freq_data = jdata[source_no][source[0]]['freq'][0]
            subtau = []
            subsf = []
            submod = {}
            submodel = {}
            subvar = {}
            max_tau_list = []
            for i in range(0, len(freq_data)):
                epoch_str = make_index_str('epoch', i)
                epoch_bin = jdata[source_no][sourcename]['epoch'][0][epoch_str]
                try:
                    max_tau = np.max(epoch_bin) - np.min(epoch_bin)
                except:
                    warnings.warn("Could not calculate maximum tau. Setting to zero.")
                    max_tau = 0
                max_tau_list.append(max_tau)
            max_tau_all = max(max_tau_list)

            # Loop through frequency bins
            for ident in range(0, len(freq_data)):
                # Creates key name for each frequency bin, for later calling
                freq_str = make_index_str('freq', ident)
                flux_str = make_index_str('flux', ident)
                epoch_str = make_index_str('epoch', ident)
                mod_str = make_index_str('mod', ident)
                numel_str = make_index_str('numel', ident)
                var_str = make_index_str('var', ident)

                # Call all data for one bin number
                freq_bin = jdata[source_no][sourcename]['freq'][0][freq_str]
                flux_bin = jdata[source_no][sourcename]['flux'][0][flux_str]
                epoch_bin = jdata[source_no][sourcename]['epoch'][0][epoch_str]
                num_epochs = jdata[source_no][sourcename]['maxNumEpochs']
                pos = jdata[source_no][sourcename]['pos']

                # Calculates the unbinned structure function and corresponding tau values
                # Assigns the same epoch to each of the previous slices made
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
                    # Create bin set for time difference, tau
                    start = 0
                    stop = max_tau
                    step = abs((start - stop)/19)

                    # Add (step) onto ending point to ensure all points are counted
                    try:
                        bin_set = list(np.arange(start, stop + step, step))
                    except ZeroDivisionError as error:
                        warnings.warn("Only 1 epoch available.")
                        bin_set = []

                    # Create empty lists for storing the binned values
                    tau_binned = [[] for _ in range(len(bin_set))]
                    sf_binned = [[] for _ in range(len(bin_set))]

                    # Bin tau and the structure function
                    [tau_binned, sf_binned] = binning2(tau, sf_list, bin_set, tau_binned, sf_binned)
                    tau_val_bin = []
                    sf_mean_bin = []

                    # Average tau and sf over each frequency bin, and append to dictionary
                    for i in bin_set:
                        halfway  = i + (1/2)*step
                        tau_val_bin.append(halfway)
                    tau_str = make_index_str('tau', ident)
                    subtau_add = {tau_str: tau_val_bin}
                    subtau.append(subtau_add)

                    for i in sf_binned:
                        if len(i) > 0:
                            av_sf = np.mean(i)
                        else:
                            av_sf = np.NaN
                        sf_mean_bin.append(av_sf)

                    sf_str = make_index_str('sf', ident)
                    subsf_add = {sf_str: sf_mean_bin}
                    subsf.append(subsf_add)
                except RuntimeWarning:
                    tau_str = make_index_str('tau', ident)
                    subtau_add = {tau_str: []}
                    subtau.append(subtau_add)
                    sf_str = make_index_str('sf', ident)
                    subsf_add = {sf_str: []}
                    subsf.append(subsf_add)

                # Try calculating the modulation index, with NaN for empty frequency bins
                # Will still return warnings
                if len(epoch_bin) > 0:
                    numel = len(np.unique(epoch_bin))
                    rms_var = np.std(flux_bin, ddof=1)
                    mean = np.mean(flux_bin)
                    mod_index = rms_var/mean
                else:
                    numel = 0
                    mod_index = np.NaN
                    warnings.warn(f"Could not calculate modulation index for frequency bin {ident}.")

                # Create dictionary structure of modulation index data, linked to
                # corresponding frequency bin
                submod_add = {mod_str: mod_index}
                submod.update(submod_add)
                submodel_add = {numel_str: numel}
                submodel.update(submodel_add)
                subvar_add = {var_str: rms_var}
                subvar.update(subvar_add)

            # Append dictionary of modulation index data to each source
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

            # Append SF and modulation index data to the empty file
            with open(kwargs['output_dir']+'/'+outname, mode='w') as f:
                outdata.append(full)
                json.dump(outdata, f)






# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculates the full structure function\
            as a function of time difference (tau), along with modulation index calculated\
            as (RMS_flux / mean_flux), which are then appended to the ouput file')
    parser.add_argument('-i', '--input_dir',
                        type=str,
                        nargs=1,
                        default=None,
                        help='location of processed data, including path')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        default=None,
                        help='location to save SF data, including path')
    args = parser.parse_args()

    main(**vars(args))
