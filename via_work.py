pip -q install 'gwosc==0.7.1'
pip -q install 'gwosc==0.7.1'
import gwosc
print(gwosc.__version__)

from gwosc.datasets import find_datasets
from gwosc import datasets

print("Current list of available catalogs")
print(find_datasets(type="catalog"))

gwtc4 = datasets.find_datasets(type='events', catalog='GWTC-4-confident')
print('GWTC-4 events:', gwtc4)
print("")

from gwosc.datasets import event_gps
gps = event_gps('GW231123')
print(gps)

from gwosc.datasets import event_at_gps
print(datasets.event_at_gps(1384782888))

from gwosc import urls
from gwosc.locate import get_event_urls
urls = get_event_urls('GW231123')
print(urls)

import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

# ! pip install -q 'gwpy==3.0.12'

import gwpy
print(gwpy.__version__)

"""`[start, end)` GPS segment to 20 seconds around this time

"""

segment = (int(gps)-10, int(gps)+10)
print(segment)

from gwpy.timeseries import TimeSeries
ldata = TimeSeries.fetch_open_data('L1', *segment, verbose = True, cache = True)
print(ldata)

plot = ldata.plot()

"""Calculating Fourier Transform of our TimeSeries"""

fft = ldata.fft()
print(fft)

from scipy.signal import get_window
window = get_window('hann', ldata.size)
lwin = ldata * window
fftamp = lwin.fft().abs()
plot = fft.abs().plot(xscale="log", yscale="log")
plot.show(warn=False)

"""Fluctuations in above plot >10 Hz seem random. This is intrinsic noise in the estimate of the spectral content of the signal from a single FFT. Reduce these fluctuations by averaging many estimates of the signal FFT; want to average are their squared moduli to find their Power Spectral Density (PSD).

Linnk to Spectral Density Guide: https://en.wikipedia.org/wiki/Spectral_density

Link to Welch's Method: https://en.wikipedia.org/wiki/Welch%27s_method
"""

asd = ldata.asd(fftlength=2, method="median")
plot = asd.plot()
plot.show(warn=False)

"""zooming in the frequency zone of interest"""

ax = plot.gca()
ax.set(xlim=(10, 1400), ylim=(1e-24, 1e-19))
plot

"""notable spikes:
*   ~80 Hz
*   ~200 Hz
*   ~400 Hz
*   ~1000 Hz

"""

ldata2 = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512, cache=True)
lasd2 = ldata2.asd(fftlength=4, method="median")
plot = lasd2.plot()
ax = plot.gca()
ax.set_xlim(10, 1400)
ax.set_ylim(1e-24, 1e-19)
plot.show(warn=False)

"""notable spikes:
*   ~80 Hz
*   ~300 Hz
*   ~500 Hz
*   ~580 Hz
*   ~610 Hz
*   ~1000 Hz
"""

#Hanford Data
hdata2 = TimeSeries.fetch_open_data('H1', int(gps)+512, int(gps)-512, cache=True)
hasd2 = hdata2.asd(fftlength=4, method="median")

#Virgo Data
#vdata2 = TimeSeries.fetch_open_data('V1', int(gps)+512, int(gps)-512, cache=True)
#vasd2 = vdata2.asd(fftlength=4, method="median")

#plot w/ standard colors
ax.plot(hasd2, label='LIGO-Hanford', color='gwpy:ligo-hanford')
#ax.plot(vasd, label='Virgo', color='gwpy:virgo')

#label and recolor Livingston Data
lline = ax.lines[0]
lline.set_color('gwpy:ligo-livingston')
lline.set_label('LIGO-Livingston')

ax.set_ylabel(r'Strain noise [$1/\sqrt{P\mathrm{Hz}}$]')
ax.legend()
plot

"""To note: Virgo data not currently available for this event."""

import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries

gps = event_gps('GW231123')
print("GW231123 GPS:", gps)

ldata = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512, cache=True)
print("GW231123 data")
print(ldata)

specgram = ldata.spectrogram2(fftlength=4, overlap=2, window='hann') ** (1/2.)
plot = specgram.plot()

ax = plot.gca()
ax.set_yscale('log')
ax.set_ylim(10, 1400)
ax.colorbar(
    clim=(1e-24, 1e-19),
    norm="log",
    label=r"Strain noise [$1/\sqrt{\mathrm{Hz}}$]",
)
plot

"""some noise noted between 30-100 Hz. Non-stationary noise."""

segment = (int(gps) - 30, int(gps) + 20)
hdata = TimeSeries.fetch_open_data('H1', *segment, verbose=True, cache=True)

hq = hdata.q_transform(frange=(30, 500))
plot = hq.plot()
plot.colorbar(label="Normalised energy")

ax = plot.gca()
ax.set_epoch(gps)
ax.set_ylim(30, 500)
ax.set_yscale("log")
plot

"""notable features at ~-14sec, ~0sec, and ~4sec"""

hq = hdata.q_transform(frange=(30, 500), qrange=(100, 110))
plot = hq.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")

hq2 = hdata.q_transform(frange=(30,500), qrange=(80, 110), outseg=(gps-3,gps+0.5))
plot = hq2.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")

hq2 = hdata.q_transform(frange=(30,500), qrange=(80, 110), outseg=(gps-5,gps+7))
plot = hq2.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")

"""NOTE: outseg intended to zoom in around merger, not currently working. (Ask Hughes)"""

hq2 = hdata.q_transform(frange=(25,100), qrange=(80, 110), outseg=(gps-0.8,gps+0.8))
plot = hq2.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")

ldata = TimeSeries.fetch_open_data('L1', *segment, verbose=True)
lq = ldata.q_transform(frange=(30, 500), qrange=(100,110))
plot = lq.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")

plot.colorbars[0].mappable.set_clim(0,20)
plot

gated_ldata = ldata.gate(tzero=0.25, tpad=0.25)

gated_lq = gated_ldata.q_transform(frange=(30, 500), qrange=(100,110))
plot = gated_lq.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")
plot.colorbars[0].mappable.set_clim(0,20)

"""# Generating Waveforms"""

import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

# ! pip install -q PyCBC==2.4.1 lalsuite==7.25

from pycbc.waveform import get_td_waveform, fd_approximants
import matplotlib.pyplot as plt

from pycbc.waveform import td_approximants

import matplotlib.pyplot as plt
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass
from pycbc.filter import get_fd_waveform

print('Time domain waveforms: ', td_approximants())
print('Frequency domain waveforms: ', fd_approximants())

merger = Merger('GW231123_135430-v2')
strain = merger.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.savefig('waveform.png')
#plt.show()

#waveform change for each mass in binary (including total mass of merger)
plt.figure(figsize=plt.figaspect(0.4))
for m in [137, 101, 223, 101]:
  hp,hc = get_td_waveform(approximant="IMRPhenomD",
                          mass1=m,
                          mass2=m,
                          delta_t=1.0/4096,
                          f_lower=101)
  plt.plot(hp.sample_times, hp, label='$M_{\odot 1,2}=%s$' % m)
plt.legend()
plt.grid

"""# Matched Filtering"""
#confirm analysis matching to GW231123 mass values
#! pip install -q 'lalsuite==7.25' 'PyCBC==2.6.0'
#import numpy
sample_rate = 1024 #samples per second
data_length = 1024 #seconds
apx = 'IMRPhenomD'

hp1, _ = get_td_waveform(approximant=apx,
                         mass1=137,
                         mass2=101,
                          delta_t=1.0/sample_rate,
                         f_lower=223)

hp1 = hp1/max(numpy.correlate(hp1, hp2, mode='full'))**0.5
plt.figure()
plt.title('The waveform hp1 for GW231123')
plt.plot(hp1.sample_times, hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')

#shift waveform to account for Guassian noise (may not be necessary for streamlined analysis)
waveform_start = numpy.random.randint(0, len[data] - len(hp1))
data[waveform_start:waveform_start+len(hp1)] += 10 * hp1.numpy()

plt.figure()
plt.plot(hp1.sample_times, data[waveform_start:waveform_start+len(hp1)])
plt.xlabel('Time(s)')
plt.ylabel('Normalized amplitude')

#plot the waveform in the noise
plt.figure()
plt.title('Signal GW231123 in the data')
plt.plot(hp1.sample_times, data[waveform_stat:waveform_start+len(hp1)])
plt.plot(hp1.sample_times, 10 * hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')
plt.show()


"""# Detection in Colored Noise"""

import pycbc.noise
import pycbc.psd

from pycbc.types import TimeSeries

#colored noise matches PSD provided (find PSD values for GW231123)
#flow = 10
#delta_f = 1.0/128 #have to check what these are
#flen = int(sample_rate / (2 * delta_f)) + 1
#psd = pycbc.psd.GW231123(flen, delta_f, flow)

#generate colored noise
delta_t = 1.0 / sample_rate
ts = pycbc.noise.noise_from_psd(data_length*sample_rate, delta_t, psd, seed = 127) #what is seed?

seg_len = int(4/ delta_t)
seg_stride = int(seg_len / 2)
estimated_psd = pycbc.psd.welch(ts, seg_len=seg_len, seg_stride=seg_stride)

plt.loglog(estimated_psd.sample_frequencies, estimated_psd, label='estimate')
plt.loglog(psd.sample_frequencies, psd, linewidth=3, label='known psd')
plt.xlim(xmin=flow, xmax=512)
plt.ylim(1e-50, 1e-43)
plt.xlabel('Frequency[Hz]')
plt.ylabel('Power spectral density')
plt.legend()
plt.grid()
plt.show()

ts[waveform_start:waveform_start+len(hp1)] += hp1.numpy() * 1E-20

cross_correlation = numpy.zeros([len(data)-len(hp1)])
hp1_numpy = hp1.numpy()
for i in range(len(data) - len(hp1_numpy)):
  cross_correlation[i] = (hp1_numpy * data[i:i+len(hp1_numpy)]).sum()

#plot cross correlated data vs. time. Superimpose location at end of signal to see where we should find peak
plt.figure()
times = numpy.arrange(len(data) - len(hp1_numpy)) / float(sample_rate)
plt.plot(times, cross_correlation)
plt.plot([waveform_start/float(sample_rate), waveform_start/float(sample_rate)], [-20,20], 'r')
plt.xlabel('Time(s)')
plt.ylabel('Cross-correlation')
plt.show()

# cross correlate the signal w entire dataset, in the time-domain

""" #Matched Filtering in Action"""
import matplotlib.pyplot as plt
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass

m = Merger('GW231123_135430-v2')
strain = merger.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.show()

""" #Filter Wraparound"""
#may not be necessary to do with our target