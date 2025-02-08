import matplotlib.pyplot as plt
import numpy as np

from funcs import binning2

def plot_comparison(output_dir, before, after):
    '''
    Plots a comparison of input data distributions before and after a cut

    Args:
        str output_dir: Location to save plot
        list before: Data values before cut
        list after: Data values after cut

    Returns:
        None.
    '''
    fig, ax = plt.subplots()

    # Create bins and plot histograms
    bin_edges = [x for x in np.linspace(min(before), max(before), endpoint=True, num=40)]
    plt.hist(before, bin_edges, facecolor='b', alpha=0.5)
    plt.hist(after, bin_edges, facecolor='r', alpha=0.5)

    # Add labels and title
    plt.xlabel('Value of modulation index at tau=50')
    plt.ylabel('Number of sources')
    plt.title('Distribution of modulation index, sampled at tau=50 before and after cut')

    # Save image and close
    plt.savefig(output_dir+'compare_modidx.pdf', format='pdf')
    plt.close()


def plot_gal(x, y, c, c_array_type, tau_val):
    '''
    Makes galactic plot with (c_array) as the data used for colour purposes

    Args:
        ndarray x: Array of x-locations of points to plot
        ndarray y: Array of y-locations of points to plot
        ndarray c: Array of color for points
        str c_type: String indicating type of data being plotted
        Numeric tau_val: Value of tau measurements were taken

    Returns:
        None
    '''
    fig, ax = plt.subplots()

    # plot dashed line for galactic plane at b=0
    plt.plot([0,0], [0,360], 'k--')

    # Scatter plot
    logc = [np.log10(ci) for ci in c]
    plt.scatter(x, y, c=logc, linewidths=0.5, cmap='gnuplot2')

    # Add colorbar and fix axis limits
    cbar = plt.colorbar()
    plt.ylim(ymax=90, ymin=-90)
    plt.xlim(xmax=360, xmin=0)

    # Set labels and savename
    plt.xlabel('l (degrees)')
    plt.ylabel('b (degrees)')
    if c_array_type == 'slope_val':
        savename = f'galplot_slope_tau_{repr(tau_val)}'
        cbar.set_label(f'log10(SF) slope at lag of {repr(tau_val)} days')
    elif c_array_type == 'fit_val':
        savename = f'galplot_fit_tau_{repr(tau_val)}'
        cbar.set_label(f'log10(SF) value at lag of {repr(tau_val)} days')
    elif c_array_type == 'mod_data':
        savename = f'galplot_mod_tau_{repr(tau_val)}'
        cbar.set_label(f'log10(modidx)')

    # Save image and close
    plt.savefig(kwargs['output_dir'] + savename, format='pdf')
    plt.close()


def plot_bar(y, var, tau):
    '''
    Creates a barplot of the distribution of var across coordinate y at lag tau

    Args:
        ndarray y: The y-location of points, given as Galactic b
        ndarray var: The variable to plot
        Numeric tau: The time lag

    Returns:
        None
    '''
    int_av = []
    int_num = []
    int_err = []

    # Create bins
    step = 5
    bin_set = list(np.arange(-90, 90 + step, step))
    b_bins = [[] for _ in range(len(bin_set))]
    int_bins = [[] for _ in range(len(bin_set))]
    b_binned, int_binned = binning2(y, var, bin_set, b_bins, int_bins)

    for i in int_binned:
        if len(i) > 0:
            int_av.append(np.nanmedian(i))
            int_num.append(len(i))
            int_err.append(statsmodels.robust.mad(i))
        else:
            int_av.append(np.NaN)
            int_num.append(0)
            int_err.append(np.NaN)

    fig, ax = plt.subplots()

    plt.bar(bin_set, int_av, step, yerr=int_err, ecolor='k', alpha=0.7)

    plt.xlim(-90,90 + step)
    plt.ylim(ymin=0)

    # Make labels and savename
    plt.xlabel('b coordinate (degrees)')
    if interest_array == fit_val:
        plt.ylabel(f'Median rate value at lag of {repr(tau)} days')
        savename = f'all_mod_fitval_5deg_median_tau_{repr(tau)}'
    elif interest_array == slope_val:
        plt.ylabel(f'Median SF slope at lag of {repr(tau)} days')
        savename = f'all_slope_5deg_median_norm_tau_{repr(tau)}'
    elif interest_array == mod_data:
        plt.ylabel(f'Median modulation index')
        savename = f'all_mod_5deg_median_norm_tau_{repr(tau)}'

    # Add bin label
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

    # Save figure and close
    plt.savefig(kwargs['output_dir'] + savename, format='pdf')
    plt.close()


def plot_coord(coords):
    '''
    Plot distribution of sources in coordinate space

    Args:
        ndarray coords: Array of point coordinates

    Returns:
        None
    '''

    bin_edges = [x for x in np.linspace(-90, 95, endpoint=True, num=38)]

    fig, ax = plt.subplots()
    plt.hist(coords, bin_edges)

    plt.xlabel('b-coordinate (degrees)')
    plt.ylabel('Number of sources')
    plt.title('Distribution of sources in b')

    plt.savefig(kwargs['output_dir'] + 'gal_dist_bcoord_5deg.pdf', format='pdf')
    plt.close()
