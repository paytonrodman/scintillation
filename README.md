# The spectral signature of interstellar scintillation

This code was originally written in the course of a summer research internship I did in the 2016-17 summer at CSIRO Astronomy and Space Science (CASS) in Perth, Australia.


## Data analysis process

The initial sample consisted of 2232 X-ray and radio quasars, sampled in frequency across a band of 2 to 12GHz. Quasars were observed using the Australian Telescope Compact Array (ATCA) between 2014 April and 2016 December, with varying periods averaging once per month. Sources were pre-selected to be free of obvious confusion, and are distributed over a wide range in Galactic coordinate space.

The full data set was obtained after calibration from the Australian Telescope Extreme Scattering Events (ATESE) pipeline, which was then read and reorganised. Frequencies below 4 and above 9GHz were infrequently observed; the considered frequency range was limited to 4.5 – 8.5GHz in order to standardise the sample. The resultant range was binned into 40 equal bins of width 100MHz and due to the possibility of edge effects, the opening and closing bins were not considered in further analysis. Information was retained for source name, sky position, epoch, and flux and frequency data, which were then outputted as a separate file for each source.

The structure function for each source is then calculated as follows. Within a single frequency bin, the magnitude of the structure function is determined by

$$ S = |F_j - F_k|^2 $$

where $F$ is the flux density measured in Jy/beam, and $j$, $k$ are two different flux elements from within that frequency bin. The time lag variable is calculated by a similar method,

$$ \tau = |t_j - t_k|, $$

for epochs $t_j$, $t_k$ in MJD format. The structure function variables are binned into 20 equal bins per frequency band, and sources with fewer than 2 epochs are automatically filtered and removed at this point. The structure function magnitude is then averaged within each bin, and the median value of $\tau$ is used in place of the mean in order to make all $\tau$ values equal across frequencies for visualisation purposes. The modulation index is calculated per frequency as

$$ m = \frac{\sqrt{\frac{1}{N} \sum_{i=1}^{N} (F_i - \overline{F})^2}}{\langle F \rangle}. $$

Data on structure function magnitude, time lags, modulation index, RMS variance of source flux density, positional information, and relevant variable sample sizes is retained and output in individual source files, which are subsequently concatenated together.

The full sample file is assessed for the purpose of flagging outliers within the modulation index data. Empty modulation index frequency bins are interpolated across by assigning them the nearest value at a lower frequency. A median window filter with window size of 15 is applied across frequencies, and the median absolute deviation (MAD) then used to determine outliers. The MAD is used in place of the standard deviation as it is less susceptible to a small number of extreme outliers. A difference between the modulation index data and the median filtered values of greater than 5 MADs is used as the cut-off. Flagged modulation index values, along with structure function information for that frequency bin, are then removed from the data set and the resulting altered data is saved to a new file. An additional log file is created detailing the specific frequencies removed by the flagging method for each source.

The structure function is then normalised by two different methods for the purpose of direct comparison between sources. The first (N1) is a simple normalisation by the mean source flux density across all epochs and frequencies:

$$ S_{N1} = \frac{|F_j - F_k|^2}{\langle F \rangle} $$

The second normalisation method (N2) converts the structure function to a rate which is preferentially weighted towards variability on shorter time-scales (on the order of 50 days), and utilises the variation in flux density for each source rather than the mean:

$$ R_{N2} = 0.02 \times \sqrt{\frac{S}{\frac{1}{N}\sum_{i=1}^{N}(F_i - \overline{F})^2}} $$

Both normalisation results are appended to the original structure function data file.

After normalisation, the sample is limited to sources with greater than 20 epochs. Within the remaining sources, modulation index values from frequency bins observed by fewer than 75% of epochs were discarded, along with the first and last modulation index bins. The structure function is sampled at ~5.7GHz, under the assumption that the structure function is well-approximated by the curve at 5.7GHz, and that either no frequencies greatly deviate from the pattern observed at 5.7GHz or that such deviations are common to all sources. Visual inspection suggests that this assumption is well founded.

The structure function is then limited to time lag values of $\le 100 days, and a linear fit is applied to the data within this range. The slope of this fit is used as a first cut to obtain the most variable sources on the right time-scales, and the particular value used as the cut-off is chosen through inspection of the distribution of slope values. In log-space, the cut-off has typically taken a value between −6 to −3 (depending on the normalisation being used) in order to return the most variable 50% of sources. This reduces the sample to on order $\sim 500$ sources. The value of the structure function along this fitted line at a time lag of $\tau = 50$ days is also determined and retained for plotting. Within this script, there are options for plotting:

 - both the reduced (with linear fit) and full structure functions,

 - comparisons in the distribution of certain values both before and after the cut is applied,

 - galactic coordinate plots of various measures of variability within a source, including slope of linear fit, value at $\tau = 50$ on linear fit, and modulation index, and

 - bar plots showing the distribution of certain variability measures as a function of b-coordinate (galactic latitude).
