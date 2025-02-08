def binning2(input1, input2, bin_set, input1_bin, input2_bin):
    '''
    Bins input1 into bin defined by the parameters of bin_set, and then averages
    over the bins.
    '''
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

def bin_frequency(freq, flux, epoch, bin_set, freq_bin, flux_bin, epoch_bin):
    '''
    Bins freq into bin defined by the parameters of bin_set, and appends
    connected values flux and epoch into separate lists with matching indices.

    Args:
        list freq: List of frequencies
        list flux: List of fluxes
        list epoch: List of epochs
        list bin_set: List of bin edges
        list of lists freq_bin: Empty placeholder lists
        list of lists flux_bin: Empty placeholder lists
        list of lists epoch_bin: Empty placeholder lists

    Returns:
        Frequencies, Fluxes, and Epochs binned according to frequency bin
    '''
    # If freq is less than a bin upper limit, place freq and the corresponding
    # flux, epoch into that bin.
    ii = 0
    jj = 0
    while ii < len(freq):
        # Reject data outside of the 4-8.5 MHz band
        if freq[ii] < 4.0 or freq[ii] > 8.5:
           ii += 1
        elif freq[ii] >= 4.0 and freq[ii] < bin_set[jj]:
           freq_bin[jj].append(freq[ii])
           flux_bin[jj].append(flux[ii])

           if not epoch[-1].isdigit():
               epoch = epoch[:-1]

           epoch_mjd = mx.DateTime.strptime(epoch, '%Y-%m-%d').mjd
           epoch_bin[jj].append(epoch_mjd)
           ii += 1
        else:
           jj += 1
    return freq_bin, flux_bin, epoch_bin


def make_index_str(name1, index):
    '''
    Creates a string for identifying the index of name1
    '''
    name1_str = name1 + repr(index)
    return name1_str


def add_ident(val_list, label):
    '''
    Adds an identifier to the beginning of each bin of input, to retain relationships
    between bins.

    Args:
        list val_list: List of values
        str label: String label to prepend

    Returns:
        list val_list: Altered list with prepended label
    '''
    m = 0
    identifier = label + repr(m)
    for i in range(len(val_list)):
       if label + repr(m) not in val_list[m]:
           val_list[m].insert(0, label + repr(m))
           m += 1
    return val_list


def append_dict(dict, sname, key, data):
    '''
    Within dictionary key 'sname', creates a dictionary with key 'type_ident',
    then appends 'add_on' under this key.

    Args:
        Dict dict:
        str sname: Source name dictionary key
        str key: String dictionary key
        list data: List of data to add to dictionary under `key`

    Returns:
        Dict dict: Altered dictionary with added keys/data
    '''
    dict[sname][key] = []
    dict[sname][key].append(data)
    return dict

def get_lb(pos):
    '''
    Converts positional data from a (RA, Dec) form to galactic coordinates (l, b)

    Args:
        dict pos: Dictionary of positional data with keys "ra" and "dec"

    Returns:
        l, b: Galactic coordinates
    '''
    ra = pos[0]['ra'].replace(':', ' ').split('.', 1)[0]
    dec = pos[0]['dec'].replace(':', ' ').split('.', 1)[0]

    coord = SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg)).galactic

    l = coord.l.degree
    b =  coord.b.degree
    return l, b
