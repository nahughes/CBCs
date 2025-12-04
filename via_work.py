pip -q install 'gwosc==0.7.1'
pip -q install 'gwosc==0.7.1'
import gwosc
print(gwosc.__version__)

#### 1.1

from gwosc.datasets import find_datasets
from gwosc import datasets

print("Current list of available catalogs")
print(find_datasets(type="catalog"))

gwtc4 = datasets.find_datasets(type='events', catalog='GWTC-4-confident')
print('GWTC-4 events:', gwtc4)
print("")

from gwosc.datasets import event_gps
gps = event_gps('GW231123_135430-v2')
print(gps)

from gwosc.datasets import event_at_gps
print(datasets.event_at_gps(1384782888)) #event gps

from gwosc.datasets import run_segment
print(run_segment('04'))

#from gwosc import urls
from gwosc.locate import get_event_urls
urls = get_event_urls('GW231123_135430-v2')
print(urls)

urls = get_event_urls('GW231123_135430-v2', duration=32, detector='H1')
print(urls)

#### 1.2

import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

import gwpy
print(gwpy.__version__)

# ! pip install -q 'gwpy==3.0.12'

segment = (int(gps)-10, int(gps)+10)
print(segment)

from gwpy.timeseries import TimeSeries
ldata = TimeSeries.fetch_open_data('L1', *segment, verbose = True, cache = True)
print(ldata)
#to read from local file use TImeseries.read method
#use dataset= function for utilizing specific dataset

plot = ldata.plot()

"""`[start, end)` GPS segment to 20 seconds around this time"""


"""Calculating Fourier Transform of our TimeSeries"""

fft = ldata.fft()
print(fft)

from scipy.signal import get_window
window = get_window('hann', ldata.size)
lwin = ldata * window
fftamp = lwin.fft().abs()
plot = fft.abs().plot(xscale="log", yscale="log")
plot.show(warn=False)

"""Fluctuations in above plot >10 Hz seem random. This is intrinsic noise 
in the estimate of the spectral content of the signal from a single FFT. 
Reduce these fluctuations by averaging many estimates of the signal FFT; 
want to average are their squared moduli to find their Power Spectral Density (PSD).

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

#Load more L1 data for more FFTs to be averaged during PSD estimation (random variations
#get averaged out)

ldata2 = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512, cache=True)
lasd2 = ldata2.asd(fftlength=4, method="median")
plot = lasd2.plot()
ax = plot.gca()
ax.set_xlim(10, 1400)
ax.set_ylim(1e-24, 1e-19)
plot.show(warn=False)

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

#### 1.3

#import warnings
#warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries

gps = event_gps('GW231123_135430-v2')
print("GW231123_135430-v2 GPS:", gps)

ldata = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512, cache=True)
print("GW231123_135430-v2 data")
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

#### 1.4
"""# Generating Waveforms"""

#import warnings
#warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

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

#waveform change for each mass in binary (including total mass of merger)
plt.figure(figsize=plt.figaspect(0.4))
for m in [137, 101, 223, 101]:
  hp,hc = get_td_waveform(approximant="IMRPhenomD",
                          mass1=m,
                          mass2=m,
                          delta_t=1.0/4096,
                          f_lower=223)
  plt.plot(hp.sample_times, hp, label='$M_{\odot 1,2}=%s$' % m)
plt.legend()
plt.grid()
plt.xlabel('Time(s)')
plt.ylabel('Strain')
plt.show()

### 2.1

"""# Intro to Matched Filtering"""
#confirm analysis matching to GW231123 mass values

#! pip install -q 'lalsuite==7.25' 'PyCBC==2.6.0'
#import numpy

apx = 'IMRPhenomD'

hp1, _ = get_td_waveform(approximant=apx,
                         mass1=137,
                         mass2=101,
                         delta_t=1.0/sample_rate,
                         f_lower=223)

hp1 = hp1/max(numpy.correlate(hp1, hp2, mode='full'))**0.5

plt.figure()
plt.title('The waveform hp1 for GW231123_135430-v2')
plt.plot(hp1.sample_times, hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')

plt.figure()
plt.title('Randomized noise')
plt.plot(hp1.sample_times, data[waveform_start:waveform_start+len(hp1)])
plt.xlabel( 'Time(s)')
plt.ylabel( 'Normalized amplitude')

plt.figure()
plt.title('Signal in the data')
plt.plot(hp1.sample_times, data[waveform_start:waveform_start+len(hp1)])
plt.plot(hp1.sample_times, 10 * hp1)
plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')
plt.show()

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

#generate a PSD for whitening the data
flow = 10.0
delta_f = 1.0 / data_length
flen = int(sample_rate / (2 * delta_f)) + 1
psd_td = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, 0)

#the PSD sampled properly for signal
delta_f = sample_rate / float(len(hp1))
flen = int(sample_rate / (2 * delta_f)) + 1
psd_hp1 = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, 0)

#0th and Nth values are zero. Set to nearby value to avoid dividing by 0
psd_td[0] = psd_td[1]
psd_td[len(psd_td) - 1] = psd_td[len(psd_td) - 2]
#do same for PSD sampled for signal
psd_hp1[0] = psd_hp1[1]
psd_hp1[len(psd_hp1) - 1] = psd_hp1[len(psd_hp1) - 2]

#convert signal and noisy data to freq domain and divide each by
#ASD=PSD ** 0.5, then convert back to time domain to whiten data
data_whitened = (ts.to_frequencyseries() / psd_td **0.5).to_timeseries()
hp1_whitened = (hp1.to_frequencyseries() / psd_hp1**0.5).to_timeseries() * 1E-21

#redo correlation in time domain with whitened data and template
cross_correlation = numpy.zeros([len(data)-len(hp1)])
hp1n = hp1_whitened.numpy()
datan = data_whitened.numpy()
for i in range(len(datan) - len(hp1n)):
  cross_correlation[i] = (hp1n * datan[i:i+len(hp1n)]).sum()

#plot cross-correlation in time domain
#superimpose location of end of signal
plt.figure()
times = numpy.arrange(len(datan) - len(hp1n)) / float(sample_rate)
plt.plot(times, cross_correlation)
plt.plot([waveform_start/float(sample_rate), waveform_start/float(sample_rate)],
            [(min(cross_correlation))*1.1, (max(cross_correlation))**1.1], 'r:')
plt.xlabel('Time(s)')
plt.ylabel('Cross-correlation')
plt.show()

#### 2.1

""""#Matched FIltering in Action""""

#same pip install as 2.1:Intro to Matched Filtering
import matplotlib.pyplot as plt
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass

merger = Merger('GW231123_135430-v2')
strain = merger.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.savefig('preconditioned waveform.png')
#plt.show()

#OR this code

m = Merger('GW231123_135430-v2')
strain = m.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.show()
plt.savefig('preconditioned waveform.png')
# cross correlate the signal w entire dataset, in the time-domain

""" #Filter Wraparound"""
#!!!may not be necessary to do with our target

#remove 2 seconds (n seconds) of data from both ends of sample time
conditioned = strain.crop(2, 2)
plt.plot(conditioned.sample_times, conditioned)
plt.xlabel('Time(s)')
plt.show()

""" #Calculate the Power Spectral Density (PSD) """
#estimate the PSD
from pycbc.psd import interpolate, inverse_spectrum_truncation

#use 4sec (n seconds) sample of time series in Welch Method (info in doc)
psd = conditioned.psd(4)

#interpolate PSD to match data and limit filter length of 1/PSD
#then directly use this PSD to filter data in a controlled manner
psd = interpolate(psd, conditioned.delta_f)

#1/PSD now act as filter with effective length of n seconds.
#Data has been highpassed above 15Hz and will have low values below this,
#so we must inform the func. to not include freq. below this freq.
psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),
                                  low_frequency_cutoff=15)

""" #Make Signal Model """

m = 223 #solar masses
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                    mass1=m,
                    mass2=m,
                    delta_t=conditioned.delta_t,
                    f_lower=223) #check on f_lower value!!!

#resize vector to match data
hp.resize(len(conditioned))

###may not be neccessary
#plot signal before and after shifting
plt.figure()
plt.title('Before shifting')
plt.plot(hp.sample_times, ho)
plt.xlabel('Time (s)')
plt.ylabel('Strain')

template = hp.cyclic_time_shift(hp.start_time)

plt.figure()
plt.title('After shifting')
plt.plot(template.sample_times, template)
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.show()

""" #Calculate Signal-to-Noise Time Series """
#need to account for length of template and 1/PSD

from pycbc.filter import matched_filter
import numpy

snr = matched_filter(template, conditioned,
                    psd=psd, low_frequency_cutoff=20)

#crop ocrrupted time(s) at beginning and end of PSD filter
snr = snr.crop( 4 + 4, 4) #seconds change based on signal ratio

plt.figure(figsize=[10,4])
plt.plot(snr.sample_times, abs(snr))
plt.ylabel('Signal-to-noise')
plt.xlabel('Time (s)')
plt.show()

peak = abs(snr).numpy().argmax()
snrp = snr[peak]
time = snr.sample_times[peak]

print("Found signal at {}s with SNR {}".format(time, abs(snrp)))

""" #Aligning and Subtracting the Proposed Signal """
from pycbc.filter import sigma

dt = time - conditioned.start_time
aligned = template.cyclic_time_shift(dt)

aligned /= signma(aligned, psd=psd, low_frequency_cutoff=20.0)

aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
aligned.start_time = conditioned.start_time

"""" #Visualize the Overlap Between the Signal and Data """
#whiten both template and data, then bandpass both between 30-300Hz

white_data = (conditioned.to_frequencyseries () / psd**0.5).to_timeseries()
white_template = (aligned.to_frequencyseries() / psd**0.5).to_timeseries()

white_data = white_data.highpass_fir(30., 512).lowpass_fir(300, 512)
white_template = white_template.highpass_fir(30, 512).lowpass_fir(300, 512)

#select time around merger
white_data = white_data.time_slice(merger.time-.2, merger.time+.1) #time around merger based on signal
white_template = white_template.time_slice(merger.time-.2, merger.time+.1)

plt.figure(figsize=[15, 3])
plt.plot(white_data.sample_times, white_data, label="Data")
plt.plot(white_template.sample_times, white_template, label="Template")
plt.legend()
plt.show()

""" #Subtracting the Signal from Data """
#now we make the chirp! :)

subtracted = conditioned - aligned

#plot the original data and the subtracted signal data

for data, title in [(conditioned, 'Original H1 Data'),
                    (subtracted, 'Signal Subtracted from H1 Data')]:

    t, f, p = data.whiten(4, 4).q_transform(.001, logfsteps=1--, qrange=(8, 8), frange=(20, 512))
    plt.figure(figsize=[15, 3])
    plt.title(title)
    plt.pcolormesh(t, f, p**0.5, vmin=1, vmax=6, shading='auto')
    plt.yscale('log')
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.xlim(merger.time - 2, merger.time + 1)
    plt.show()


####2.3

""" Read and Precondition GW Strain Data """
#in this unit we will precondition as before, calculate PSD, and re-weight SNR
#time-series of different detector strains

#useful for cross-examining strains from multiple detectors

m = Merger("GW231123_135430")
ifos = ['H1', 'L1']
data = {}
psd = {}

for ifo in ifos:
    ts = m.strain(ifo).highpass_fir(15, 512)
    data[ifo] = resample_to_delta_t(ts, 1.0/2048).crop(2, 2)

    p = data[ifo].psd(2)
    p = interpolate(p, data[ifo].delta_f)
    p = inverse_spectrum_truncation(p, int(2 * data[ifo].sample_rate), low_frequency_cutoff=15.0)
    psd[ifo] = p 
    plt.plot(psd[ifo].sample_frequencies, psd[ifo], label=ifo)

  plt.yscale('log')
  plt.xscale('log')
  plt.ylim(1e-47, 1e-43)
  plt.xlim(20, 1024)
  plt.ylabel('$Strain^2 / Hz$')
  plt.xlabel('Frequency (Hz)')
  plt.grid()
  plt.legend()
  plt.show()

""" #Generate waveform and Calculate Sig-to-noise Time Series """
#calc component mass of each BH in frame
cmass = (m.median1d("mass1")+m.median1d("mass2")) / 2
cmass *= (1 + m.median1d("redshift"))

hp, _ = get_fd_waveform(approximant="IMRPhenomD",
                        mass1=cmass, mass2=cmass,
                        f_lower=20.0, delta_f=data[ifo].delta_f)
hp.resize(len(psd[ifo]))

snr={}
for ifo in ifos:
    snr[ifo] = matched_filter(hp, data[ifo], psd=psd[ifo], low_frequency_cutoff=20)
    snr[ifo] = snr[ifo].crop(5, 4) #crop based on signal

#various sizes
for w, title in [(8, 'Wide View'), (.15, 'Close to GW231123')]:
    plt.figure(figsize=[14, 4])
    for ifo in ifos:
        plt.plot(snr[ifo].sample_times, abs(snr[ifo]), label=ifo)
    
    plt.legend()
    plt.title(title)
    plt.grid()
    plt.xlim(m.time - w, m.time + w)
    plt.ylim(0, 15)
    plt.xlabel('Time (s)')
    plt.ylabel('Signal-to-noise (SNR)')
    plt.show()

#calculate SNR of p bins to the total SNR:
# x^2 = p sigma i=0 (p_i - p/p)^2
#will have 2p-2 degrees of freedom (each SNR is complex, both possible orthogonal phases signal
#could have contributions from)
#constraint due to sum of each bin must add up to total SNR
#normalize this by dividing by number of degrees of freedom, producing X_r^2

#if having issues with following code, replace import w/ from_pycbc_chisq import power_chisq
from pycbc.vetoes import power_chisq

chisq = {}
for ifo in ifos:
    #nbins arbitrary amount
    nbins = 30
    chisq[ifo] = power_chisq(hp, data[ifo], nbins, psd[ifo], low_frequency_cutoff=20.0)
    chisq[ifo] = chisq[ifo].crop(5, 4)

    dof = nbins * 2 - 2
    chisq[ifo] /= dof

    plt.figure(figsize=[14, 4])
    for ifo in ifos:
        plt.plot(chisq[ifo].sample_times, chisq[ifo], label=ifo)

    plt.legend()
    plt.grid()
    plt.xlim(m.time -0.15, m.time + 0.15)
    plt.ylim(0, 5)
    plt.ylabel('$chi^2_r$')
    plt.show()

""" #Reweight SNR to down-weight times not in signal """
#combine SNR time series and our X_r^2 time series
#two papers explaining process in Doc

from pycbc.events.ranking import newsnr

#rho-hat symbol above (in paper) named "newsnr" here
nsnr = {ifo:newsnr(abs(snr[ifo])), chisq[ifo] for ifo in ifos}

for w, title in [(8, 'Wide View'), (.15, 'Close to GW231123')]:
    plt.figure(figsize[14, 4])
    for ifo in ifos:
        plt.plot(snr[ifo].sample_times, nsnr[ifo], label=ifo)
    
    plt.legend()
    plt.title(title)
    plt.grid()
    plt.xlim(m.time - w, m.time + w)
    plt.ylim(0, 15)
    plt.xlabel('Time (s)')
    plt.ylabel('Re-weighted SNR')
    plt.show()


""" #Calculating Background and Significance """
#calculate light-travel-time between detectors
from pycbc.detector import Detector

#tof = "time of flight"
d = Detector("H1")
tof = {}
tof['L1'] = d.light_travel_time_to_detector(Detector("L1"))

#record time of peak in LIGO Observatories
ptime = {}

#shade region around each LIGO peak that could have peak in Hanford if from astrophysical source
plt.figure(figsize=[14, 4])
for ifo in ifos:
    if ifo != 'H1':
        ptime[ifo] = snr[ifo].sample_times[nsnr[ifo].argmax()]
        plt.axvspan(ptime[ifo] - tof[ifo], ptime[ifo] + tof[ifo], alpha=0.2, lw=10)
    plt.plot(snr[ifo].sample_times, nsnr[ifo], label=ifo)

#ISSUE: need to calculate span of time that Hanford peak could happen compared to
#both other detectors (Virgo not available for GW231123)

# start = ptime['H1] - tof['H1]
# end = ptime['L1] + tof ['L1]

# convert times to indices along how large region is
# window_size = int((end-start) * snr['H1].sample_rate)
# sidx = int((start - snr['H1].start_time) * snr['L1'].sample_rate)
# eidx = sidx + window_size

#Calculate "on-source" peak re-weighted (newsnr) statistic value
onsource = nsnr['H1'][sidx:eidx].max()

plt.legend()
plt.grid()
plt.xlim(m.time - .08, m.time + .08)
plt.ylim(0, 15)
plt.xlabel('Time (s)')
plt.ylabel('Re-weighted SNR')
plt.show()

print('Hanford peak has a re-weighted SNR value of {}'.format(onsource))

#walk through data in chunks to calc. peak stat. value in each
peaks = []
i = 0
while i + window_size < len(nsnr['H1']):
    p = nsnr['H1'][i:i+window_size].max()
    peaks.append(p)
    i += window_size
    if abs(i - sidx) < window_size:
        i += window_size * 2
peaks = numpy.array(peaks)

#make mapping btwn stat. value and p-value (number of observed background samples w/ value
# >= the onsource divided by number of samples) using our background samples
pcurve = numpy.arange(1, len(peaks)+1)[::-1] / float(len(peaks))
peaks.sort()

pvalue = (peaks > onsource).sum() / float(len(peaks))

plt.figure(figsize=[10, 7])
plt.scatter(peaks, pcurve, label='Off-source (Noise Background)', color='black')

plt.axvline(onsource, label='On-source', color='red')
plt.axhline(pvalue, color='red')

plt.legend()
plt.yscale('log')
plt.grid()
plt.ylim(1e-5, 1e0)
plt.ylabel('p-value')
plt.xlabel('Re-weighted SNR')

plt.xlim(1, 6)
plt.show()
print("The p-value associated w/ GW231123 peak is {}".format(pvalue))

""" #Introduction to Parameter Estimation """

#guide use example model, we will attempt w/ GW231123 data
