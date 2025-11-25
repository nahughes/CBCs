import numpy
# The first import of matplotlib can take some time (especially on cloud platforms). This is normal.
import matplotlib.pyplot as plt

import gwosc
# print("GWOSC Version")
# print(gwosc.__version__)
# print("")

from gwosc.datasets import find_datasets
from gwosc import datasets
# print("Current list of available catalogs")
# print(find_datasets(type="catalog"))
gwtc4 = datasets.find_datasets(type='events', catalog='GWTC-4.0')
# print('GWTC-4 events:', gwtc4)
# print("")

from gwosc.datasets import event_gps
gps = event_gps('GW231123_135430-v2')
print("Event GPS")
print(gps)
print("")

from gwosc.locate import get_event_urls
urls = get_event_urls('GW231123_135430-v2')
# print(urls)

import gwpy
# print("GWPY Version")
# print(gwpy.__version__)
# print("")

segment = (int(gps)-512, int(gps)+512)
print("Timeframe")
print(segment)
print("")

from gwpy.timeseries import TimeSeries
ldata = TimeSeries.fetch_open_data('L1', *segment, verbose = True, cache = True)
print("Livingston Strain")
print(ldata)
print("")

# specify the sample rate.
# LIGO raw data is sampled at 16384 Hz (=2^14 samples/second).
# It captures signal frequency content up to f_Nyquist = 8192 Hz.
# Here, we will make the computation faster by sampling at a lower rate.
sample_rate = 4096 # samples per second
data_length = len(ldata) # seconds
data = ldata

# ldata.plot()
# plt.savefig('ligo-strain-test.png')
# print("plot done!")

# plt.plot(times, data)
data.plot()
plt.xlabel('Time(s)')
plt.savefig('waveform-rawdata.png')
print("Strain data plot done!")

# data = numpy.random.normal(size=[sample_rate * data_length])
times = numpy.arange(len(data)) / float(sample_rate)

from pycbc.waveform.waveform import get_td_waveform

# the "approximant" (jargon for parameterized waveform family).
# IMRPhenomD(a phenomenological Inspiral–Merger–Ringdown wafeform model) is defined in the frequency domain, but we'll get it in the time domain (td).
# It runs fast, but it doesn't include effects such as non-aligned component spin, or higher order modes.
apx = 'IMRPhenomD'

hp1, _ = get_td_waveform(approximant=apx,
                         mass1=137,
                         mass2=101,
                          delta_t=1.0/sample_rate,
                         f_lower=2)

# hp1 = hp1 / max(numpy.correlate(hp1, hp1, mode='full'))**0.5

# note that in this figure, the waveform amplitude is of order 1.
# The duration (for frequency above f_lower=25 Hz) is only 3 or 4 seconds long.
# The waveform is "tapered": slowly ramped up from zero to full strength, over the first second or so.
# It is zero-padded at earlier times.
plt.figure()
plt.title("The waveform hp1")
plt.plot(hp1.sample_times, hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')
plt.savefig('waveform_approx.png')
print("Waveform plot done!")


plt.figure()
plt.title("Signal in the data")
plt.plot(times, data)
plt.plot(hp1.sample_times, 10 * hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')
plt.savefig('waveform-overlay.png')
print("Overlay plot done!")

## Cross Correlation is currently taking a very long time to complete.
'''
cross_correlation = numpy.zeros([len(data)-len(hp1)])
hp1_numpy = hp1.numpy()
for i in range(len(data) - len(hp1_numpy)):
    cross_correlation[i] = (hp1_numpy * data[i:i+len(hp1_numpy)]).sum()

# plot the cross-correlated data vs time. Superimpose the location of the end of the signal;
# this is where we should find a peak in the cross-correlation.
plt.figure()
times = numpy.arange(len(data) - len(hp1_numpy)) / float(sample_rate)
plt.plot(times, cross_correlation)
plt.plot([waveform_start/float(sample_rate), waveform_start/float(sample_rate)], [-10,10],'r:')
plt.xlabel('Time (s)')
plt.ylabel('Cross-correlation')
plt.savefig('cross_correl.png')
'''