#!/usr/bin/python
#
# make_plots.py
#
# A program to determine which sources are likely experiencing high instrinic
# variability based on their Structure function, and then plot them
#
# Usage: python make_plots.py -i [input_dir] -o [output_dir]
#
# Output:
#    - a .json file with the Structure Function and modulation index
#

# Python standard modules
import argparse
import sys
import warnings

# Other Python modules
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
import json
import numpy as np
import statsmodels.robust

from funcs import binning2, make_index_str, get_xy
import plot_funcs


def main(**kwargs):
    # Read SF data
    data = []
    with open(kwargs['input_file']) as f:
        for line in f:
            data.append(json.loads(line))
    data = data[0][0] #get through 2 levels of list wrapping

    # Initialise empty lists/dicts for later
    x_list = []
    y_list = []
    mod_data = []
    outfile = []
    fit_val = []
    slope_val = []
    slope_all = []
    mod_all = []
    full = {}
    passcounter = 0

    #create set of frequency bins
    start = 4.5
    stop = 8.5
    step = (stop - start)/39 #must be same as data_extract.py -> check this
    bin_set = list(np.arange(start, stop + step, step)) #must add 1*(step) to desired end-point.

    for source_no in range(len(data)):
        source_name = list(data[source_no][0].keys())[0]

        ddict = data[source_no][0][source_name]
        if ddict['num_epochs'] > 20:
            print(f'Source number [{source_no+1}] out of [{len(data)}]')

            # Read in positional data and convert RA, Dec from (hh:mm:ss) and degrees
            # into galactic coordinates (l,b)
            l, b = get_lb(data_type['pos'])

            freq_list = []
            mod_list = []
            for j in range(len(ddict['tau'][0])):
                num_epoch_meas = ddict['mod_numel'][0][make_index_str('numel', j)]

                # Make a list of all modulation indices, for comparison to the result from the cut
                if j == 0 or j == len(ddict['tau'][0]) - 1: # discard first and last bins
                    mod_el = np.NaN
                else:
                    if num_epoch_meas >= (0.75*ddict['num_epochs']): # if observed by >75% of epochs, retain
                        mod_el = ddict['mod'][0][make_index_str('mod', j)]
                    else: # else, discard
                        mod_el = np.NaN
                mod_list.append(mod_el)

            # Loop through frequencies, not including the first and last bins
            # Take frequency bin at ~5.7GHz as approximation for all frequencies ('midway' frequency)
            midway = 5
            for j in range(1, len(ddict['tau'][0]) - 1):
                num_epoch_meas = ddict['mod_numel'][0][make_index_str('numel', j)]
                tau_bin = ddict['tau'][0][j][make_index_str('tau', j)]
                sf_bin = ddict['sf'][0][j][make_index_str('sf', j)]
                rate_bin = ddict['rate'][0][make_index_str('rate', j)]
                normsf_bin = ddict['sf_norm'][0][make_index_str('sfnorm', j)]

                if num_epoch_meas >= (0.75*ddict['num_epochs']): # if observed by >75% of epochs, keep
                    mod_val = ddict['mod'][0][make_index_str('mod', j)]
                else: # else, discard
                    mod_val = np.NaN

                # Get frequency for this index (nominally 5.7GHz)
                try:
                    midfreq = bin_set[midway - 1] # (-1) as midway relates to length and is not indexed from 0
                except:
                    midfreq = np.NaN

                # If frequency is ~5.7GHz, find all tau values that are less
                # than 100, and corresponding SF values
                if j == midway:
                    lim = 100 # limit in tau that we fit to
                    indices = [z for z,val in enumerate(tau_bin) if val <= lim] # indices of all tau points that are <= lim
                    tau_all = []
                    sf_all = []

                    for i in indices:
                        tau_all.append(tau_bin[i])
                        sf_all.append(sf_bin[i])

                    # Make a linear fit to tau, and plot
                    coeff = np.polyfit(tau_all, sf_all, 1) # non-normalized structure function
                    fit = np.poly1d(coeff)
                    plt.plot(tau_all, fit(tau_all))
                    slope = coeff[0]
                    tau_sample = 50 # tau value to sample at
                    testing_val = (slope*tau_sample) + coeff[1] # value at tau = tau_sample

                    if slope < 0: # negative values fail when taking the log
                        slope = 10**(-16)

                    cutoff = -6 # cutoff for selecting approx top 50% of sources
                    mod_all.append(mod_val)

                    # Use log10(slope) of linear fit as an indication for
                    # modulation index (log10 as slope becomes saturated otherwise)
                    if np.log10(slope) >= cutoff:
                        passcounter += 1 # counts the number of sources which pass the cut
                        fit_val.append(testing_val)
                        slope_val.append(slope)
                        mod_data.append(mod_val)
                        b_list.append(b)
                        l_list.append(l)
                        full[source_name] = mod_list

            # Update plot qualities, save and then close
            plt.xlabel('Binned time difference (tau, MJD)')
            plt.ylabel('Binned Structure Function (Jy/beam, squared)')
            plt.title(f'Source name: {source_name}, Number of epochs: {repr(ddict['num_epochs']}')
            ax = plt.gca()
            ax.set_ylim(ymin=0)
            plt.savefig('sf_' + source_name, format='pdf', bbox_inches='tight')
            plt.close()

    print(f'{passcounter} sources passed the cut')


    plot_funcs.plot_comparison(mod_all, mod_data)
    plot_funcs.plot_gal(l_list, b_list, slope_val, 'slope_val', tau_sample)
    plot_funcs.plot_gal(l_list, b_list, fit_val, 'fit_val', tau_sample)
    plot_funcs.plot_gal(l_list, b_list, mod_data, 'mod_data', tau_sample)
    plot_funcs.plot_bar(b_list, fit_val, tau_sample)
    plot_funcs.plot_bar(b_list, slope_val, tau_sample)
    plot_funcs.plot_bar(b_list, mod_data, tau_sample)
    plot_funcs.plot_coord(b_list)

    log = [np.log10(x) for x in slope_val]
    bin_edges = [x for x in np.linspace(min(log), max(log), endpoint=True, num=40)]
    plt.hist(log, bin_edges, facecolor='b', alpha=0.5)
    plt.xlabel('Distribution of log10(slope) of rate function')
    plt.ylabel('Number of sources')
    plt.show()

    #save new data structure
    with open(kwargs['output_dir']+'/variable_sources.json', mode='w') as empty:
        outfile.append(full)
        json.dump(outfile, empty)


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine type of variability\
    based on Structure Function slope at a time lag of tau=50')
    parser.add_argument('-i', '--input_file',
                        type=str,
                        nargs=1,
                        default=None,
                        help='input structure function file, including path')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        nargs=1,
                        default=None,
                        help='location to save output files, including path')
    args = parser.parse_args()

    main(**vars(args))
