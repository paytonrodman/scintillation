#!/usr/bin/python
#
# normalize_sf.py
#
# A program to normalize the cleaned Structure Function for each source.
#
# Usage: python normalize_sf.py -i [input_file] -s [source_dir] -o [output_dir]
#
# Output:
#    - a .json file containing the normalized Structure Functions
#

# Python standard modules
import argparse
import sys
import glob

# Other Python modules
import json
import numpy as np

def main(**kwargs):
    # Get all Structure Function data
    sf_all = []
    with open(kwargs['input_file']) as f:
        for line in f:
            sf_all.append(json.loads(line))

    for source_no in range(len(sf_all[0])):
        subrate = {}
        subsfnorm = {}
        sourcename = sf_all[0][source_no][0].keys()
        sourcename = sourcename[0]
        print(f'Source number [{source_no+1}] out of [{len(sf_all[0])}] in sf_all')

        sf = sf_all[0][source_no][0][sourcename]['sf']

        # Get flux data for source from original json
        f_data = json.load(open(kwargs['source_dir'] + sourcename + '_SF.json'))
        flux_all = []
        for freq_bin in f_data[0][sourcename]['flux'][0]:
            flux_all.append(f_data[0][sourcename]['flux'][0][freq_bin])

        av_flux = np.mean(flux_all[0])
        sf_all[0][source_no][0][sourcename]['av_flux'] = av_sourceflux
        rms = sf_all[0][source_no][0][sourcename]['rms_var']

        # Normalise SF by two different methods (rate and sf_norm)
        for freq_num in range(len(sf[0])):
            sf_bin = sf[0][freq_num]['sf_'+repr(freq_num)]
            rms_bin = rms[0]['var_'+repr(freq_num)]

            try:
                rate_bin_norm = [(0.02*(np.sqrt(x/(rms_bin**2)))) for x in sf_bin]
            except ZeroDivisionError as error:
                rate_bin_norm = [np.NaN for x in sf_bin]

            sf_bin_norm = [(np.sqrt(x)/av_sourceflux) for x in sf_bin]

            subrate_add = {'rate_'+repr(freq_num): rate_bin_norm}
            subrate.update(subrate_add)
            subsfnorm_add = {'sfnorm_'+repr(freq_num): sf_bin_norm}
            subsfnorm.update(subsfnorm_add)

        # Append normalisation to SF data
        sf_all[0][source_no][0][sourcename]['rate'] = []
        sf_all[0][source_no][0][sourcename]['rate'].append(subrate)
        sf_all[0][source_no][0][sourcename]['sf_norm'] = []
        sf_all[0][source_no][0][sourcename]['sf_norm'].append(subsfnorm)

    # Output normalisation in new file
    with open(kwargs['output_dir']+'SF_normalized.', mode='w') as f:
        f.write(json.dumps(sfdata))


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_dir',
                        type=argparse.FileType('r'),
                        nargs=1,
                        default=None,
                        help='Cleaned structure function data file, including path.')
    parser.add_argument('-s', '--source_dir',
                        type=str,
                        default=None,
                        help='Directory containing source files')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        default=None,
                        help='Directory to save output')
    args = parser.parse_args()

    main(**vars(args))
