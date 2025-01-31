#!/usr/bin/python

# import all necessary packages
import matplotlib
from matplotlib import pyplot
import astropy
from astropy import units
from astropy.coordinates import SkyCoord
import json
import numpy as np
import sys
import statsmodels.robust

''' Bins input1 into bin defined by the parameters of bin_set, and then averages
over the bins'''
def binning(input1, input2, bin_set, input1_bin, input2_bin):
    #if bcoord is less than first bin limit, place in first bin otherwise, test
    #if bcoord is between two bins, and place in corresponding bin
    l = 0
    for i in input1:
        if input1[l] <= bin_set[0]:
            input1_bin[0].append(input1[l])
            input2_bin[0].append(input2[l])
        else:
            for k in range(1, len(bin_set)):
                if input1[l] <= bin_set[k] and input1[l] > bin_set[k - 1]:
                    input1_bin[k].append(input1[l])
                    input2_bin[k].append(input2[l])
        l += 1
    return input1_bin, input2_bin

''' Creates a string for identifying the bin number (index) of name1'''
def bin_str(name1, index):
    name1_str = name1 + repr(index)
    return name1_str

'''Converts positional data from a RA+Dec form to galactic coordinates (l,b)'''
def getxy(pos):
    ra = pos[0]['ra']
    dec = pos[0]['dec']
    ra =  ra.replace(':', ' ')
    dec =  dec.replace(':', ' ')
    sep  = '.'
    ra = ra.split(sep, 1)[0]
    dec = dec.split(sep, 1)[0]
    pos_deg = ra + ' ' + dec
    coord = SkyCoord(pos_deg, unit=(units.hourangle, units.deg))
    coord = coord.galactic
    y =  coord.b.degree
    x = coord.l.degree
    return x,y

'''Plots before and after cut distributions together for input data'''
def comparison(input_before, input_after):
    bin_edges = [x for x in np.linspace(min(input_before), max(input_before), endpoint=True, num=40)]
    pyplot.hist(input_before, bin_edges, facecolor='b', alpha=0.5)
    pyplot.hist(input_after, bin_edges, facecolor='r', alpha=0.5)
    pyplot.xlabel('Value of modulation index at tau=50')
    pyplot.ylabel('Number of sources')
    pyplot.title('Distribution of modulation index, sampled at tau=50 before and after cut')
    outfilename = 'compare_modidx'

    pyplot.savefig(outfilename, format='pdf')
    #pyplot.show()
    pyplot.close()

'''Makes galactic plot with (c_array) as the data used for colour purposes'''
def galplot(x_array, y_array, c_array, c_array_type, tau):
    #plot dashed line for galactic plane at b=0
    yval = [0,0]
    xval = [0,360]
    pyplot.plot(xval, yval, 'k--')
    c_array = [np.log10(x) for x in c_array]
    pyplot.scatter(x_array, y_array, c=c_array, linewidths=0.5, cmap='gnuplot2')
    pyplot.xlabel('l (degrees)')
    pyplot.ylabel('b (degrees)')
    cbar = pyplot.colorbar()
    pyplot.ylim(ymax=90, ymin=-90)
    pyplot.xlim(xmax=360, xmin=0)
    if c_array_type == 'slope_val':
        savename = 'galplot_slope_tau' + repr(tau)
        cbar.set_label('log10(SF) slope at lag of ' + repr(tau) + ' days')
    elif c_array_type == 'fit_val':
        savename = 'galplot_fit_tau' + repr(tau)
        cbar.set_label('log10(SF) value at lag of ' + repr(tau) + ' days')
    elif c_array_type == 'mod_data':
        savename = 'galplot_mod_tau' + repr(tau)
        cbar.set_label('log10(modidx)')
    pyplot.savefig(savename, format='pdf')
    #pyplot.show()
    pyplot.close()

'''Creates a barplot of the distribution of (interest_array) across coordinates
in (y_array), and at lag (tau)'''
def barplt(y_array, interest_array, tau):
    int_av = []
    int_num = []
    int_err = []
    step = 5
    bin_set = list(np.arange(-90, 90 + step, step))
    b_bins = [[] for _ in range(len(bin_set))]
    int_bins = [[] for _ in range(len(bin_set))]
    b_binned, int_binned = binning(y_array, interest_array, bin_set, b_bins, int_bins)
    for i in int_binned:
        if len(i) > 0:
            int_av.append(np.nanmedian(i))
            int_num.append(len(i))
            #int_err.append(np.nanstd(i, ddof=1))
            int_err.append(statsmodels.robust.mad(i))
        else:
            int_av.append(np.NaN)
            int_num.append(0)
            int_err.append(np.NaN)
    bar_width = step
    pyplot.bar(bin_set, int_av, bar_width, yerr=int_err, ecolor='k', alpha=0.7)
    #pyplot.bar(bin_set, int_av, bar_width, alpha=0.7)
    ax = pyplot.gca()
    pyplot.xlim(-90,90 + step)
    pyplot.xlabel('b coordinate (degrees)')
    if interest_array == fit_val:
        pyplot.ylabel('Median rate value at lag of ' + repr(tau) + ' days')
        savename = 'all_mad_fitval_5deg_median_tau' + repr(tau)
    elif interest_array == slope_val:
        pyplot.ylabel('Median SF slope at lag of ' + repr(tau) + ' days')
        savename = 'all_slope_5deg_median_norm_tau' + repr(tau)
    elif interest_array == mod_data:
        pyplot.ylabel('Median modulation index')
        savename = 'all_mod_5deg_median_norm_tau' + repr(tau)
    pyplot.ylim(ymin=0)
    rects = ax.patches
    labels  = [repr(i) for i in int_num if not np.isnan(i)]
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom
    for rect,label in zip(rects,labels):
        height = rect.get_height()
        if np.isnan(height):
            height = 0
        label_pos = height + (y_height*0.01)
        ax.text(rect.get_x() + rect.get_width()/2, label_pos, label, ha='center', va='bottom')
    pyplot.savefig(savename, format='pdf')
    #pyplot.show()
    pyplot.close()

'''Plot distribution of sources in coordinate space'''
def coord_dist(coord_array):
    bin_edges = [x for x in np.linspace(-90, 95, endpoint=True, num=38)]
    pyplot.hist(coord_array, bin_edges)
    pyplot.xlabel('b-coordinate (degrees)')
    pyplot.ylabel('Number of sources')
    pyplot.title('Distribution of sources in b')
    pyplot.savefig('gal_dist_bcoord_5deg', format='pdf')
    #pyplot.show()
    pyplot.close()




'''Reads pre-calculated structure function and reduces this to the first (1/3)
fraction of tau values. A linear fit is made to tau values between 0 and 100, and
the value of the slope of this fit at tau=50 is used as a measure of the likelihood
to be an intrinsic variable or otherwise. The reduced structure function and linear
fit are plotted and saved for all sources. If the slope is above some arbitrary
cutoff, then the source is added to a galactic position plot, which is then saved.'''
infilename = sys.argv[1] #infilename should be the outfile from 'calculate_sf.py' - all files concatenated
outfilename = sys.argv[2] #file of modulation indices for sources above sf cutoff

data = []
with open(infilename) as f:
    for line in f:
        data.append(json.loads(line))
data = data[0][0] #get through 2 levels of list wrapping
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
    sourcename = data[source_no][0].keys()
    sourcename = sourcename[0]
    if sourcename in data[source_no][0].keys(): #optional source selection method
        data_type = data[source_no][0][sourcename] #gives dictionary to ['tau', 'sf', 'mod', 'rate', 'sf_norm', 'num_epochs', 'mod_numel']
        tau = data_type['tau']
        rate = data_type['rate']
        sf = data_type['sf']
        normsf = data_type['sf_norm']
        mod = data_type['mod']
        num_epochs = data_type['num_epochs']
        modnumel = data_type['mod_numel']
        if num_epochs > 20:
            print 'Source number [', source_no + 1, '] out of [', len(data), ']'
            #read in positional data and convert RA and Dec from (hh:mm:ss) and degrees
            # into galactic coordinates (l,b)
            pos = data_type['pos']
            x, y = getxy(pos)
            freq_list = []
            mod_list = []
            for j in range(len(tau[0])):
                mod_key = bin_str('mod', j)
                num_key = bin_str('numel', j)
                num_epoch_meas = modnumel[0][num_key]
                #make a list of all modulation indices, for comparison to the result from the cut
                if j == 0 or j == len(tau[0]) - 1: #discard first and last bins
                    mod_el = np.NaN
                else:
                    if num_epoch_meas >= (0.75*num_epochs): #if observed by >75% of epochs, retain
                        mod_el = mod[0][mod_key]
                    else: #else, discard
                        mod_el = np.NaN
                mod_list.append(mod_el)
            #loop through frequencies, not including the first and last bins
            #take frequency bin at ~5.7GHz as approximation for all frequencies ('midway' frequency)
            midway = 5
            for j in range(1, len(tau[0]) - 1):
                #keys for accessing each frequency bin
                tau_key = bin_str('tau', j)
                sf_key = bin_str('sf', j)
                rate_key = bin_str('rate', j)
                normsf_key = bin_str('sfnorm', j)
                mod_key = bin_str('mod', j)
                num_key = bin_str('numel', j)

                #accessing the bins
                num_epoch_meas = modnumel[0][num_key]
                tau_bin = tau[0][j][tau_key]
                sf_bin = sf[0][j][sf_key]
                rate_bin = rate[0][rate_key]
                normsf_bin = normsf[0][normsf_key]
                if num_epoch_meas >= (0.75*num_epochs): #if observed by >75% of epochs, keep
                    mod_val = mod[0][mod_key]
                else: #else, discard
                    mod_val = np.NaN
                try: #get frequency for this index (nominally 5.7GHz)
                    midfreq = bin_set[midway - 1] #(-1) as midway relates to length and is not indexed from 0
                except:
                    midfreq = np.NaN
                if j == midway: #if frequency is ~5.7GHz
                    #find all tau values that are less than 100, and corresponding
                    # sf values
                    lim = 100 #limit in tau that we fit to
                    indices = [z for z,val in enumerate(tau_bin) if val <= lim] #indices of all tau points that are <= lim
                    tau2fit = []
                    sf2fit = []
                    normsf2fit = []
                    rate2fit = []
                    for i in indices:
                        tau2fit.append(tau_bin[i])
                        sf2fit.append(sf_bin[i])
                        normsf2fit.append(normsf_bin[i])
                        rate2fit.append(rate_bin[i])
                    #make a linear fit to tau2fit, and plot
                    try:
                        coeff = np.polyfit(tau2fit, sf2fit, 1) #unnormalized structure function
                        #coeff = np.polyfit(tau2fit, normsf2fit, 1) #N2 normalized
                        #coeff = np.polyfit(tau2fit, rate2fit, 1) #N1 normalized
                    except: #if the bin at ~5.7GHz is empty, then try the next bin (doesn't ever occur)
                        midway += 1
                        print 'sampling at midway = ', midway
                        continue

                    fit = np.poly1d(coeff)
                    pyplot.plot(tau2fit, fit(tau2fit))
                    slope = coeff[0]
                    val_at = 50 #tau value to sample at
                    testing_val = (slope*val_at) + coeff[1] #value at tau = val_at

                    if slope < 0: #negative values fail when taking the log
                        slope = 10**(-16)
                    cutoff = -6 #cutoff for selecting approx top 50% of sources
                    mod_all.append(mod_val)
                    if np.log10(slope) >= cutoff: #use log10(slope) of linear fit as an indication for modulation index (log10 as slope becomes saturated otherwise)
                        passcounter += 1 #counts the number of sources which pass the cut
                        fit_val.append(testing_val)
                        slope_val.append(slope)
                        mod_data.append(mod_val)
                        y_list.append(y)
                        x_list.append(x)
                        full[sourcename] = mod_list
            #update plot qualities, save and then close
            pyplot.xlabel('Binned time difference (tau, MJD)')
            pyplot.ylabel('Binned Structure Function (Jy/beam, squared)')
            string1 = 'Source name: ' + sourcename + ', Number of epochs: ' + repr(num_epochs)
            pyplot.title(string1)
            ax = pyplot.gca()
            ax.set_ylim(ymin=0)
            #pyplot.savefig('sf_' + sourcename, format='pdf', bbox_inches='tight')
            pyplot.close()
print passcounter, 'sources passed the cut'


comparison(mod_all, mod_data)
galplot(x_list, y_list, slope_val, 'slope_val', val_at)
galplot(x_list, y_list, fit_val, 'fit_val', val_at)
galplot(x_list, y_list, mod_data, 'mod_data', val_at)
barplt(y_list, fit_val, val_at)
barplt(y_list, slope_val, val_at)
barplt(y_list, mod_data, val_at)
coord_dist(y_list)

log = [np.log10(x) for x in slope_val]
bin_edges = [x for x in np.linspace(min(log), max(log), endpoint=True, num=40)]
pyplot.hist(log, bin_edges, facecolor='b', alpha=0.5)
pyplot.xlabel('Distribution of log10(slope) of rate function')
pyplot.ylabel('Number of sources')
pyplot.show()

#save new data structure
with open(outfilename, mode='w') as empty:
    outfile.append(full)
    json.dump(outfile, empty)
