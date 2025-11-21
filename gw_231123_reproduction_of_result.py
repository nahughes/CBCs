import matplotlib.pyplot as plt
# import numpy as np
# import subprocess
# import shlex

import gwosc
print(gwosc.__version__)
print("")

from gwosc.datasets import find_datasets
from gwosc import datasets
# print("Current list of available catalogs")
# print(find_datasets(type="catalog"))
gwtc4 = datasets.find_datasets(type='events', catalog='GWTC-4.0')
print('GWTC-4 events:', gwtc4)
print("")

from gwosc.datasets import event_gps
gps = event_gps('GW231123_135430-v2')
print(gps)

from gwosc.locate import get_event_urls
urls = get_event_urls('GW231123_135430-v2')
print(urls)

import gwpy
print(gwpy.__version__)

segment = (int(gps)-512, int(gps)+512)
print(segment)

from gwpy.timeseries import TimeSeries
ldata = TimeSeries.fetch_open_data('L1', *segment, verbose = True, cache = True)
print(ldata)
print("")
plot = ldata.plot()
plt.savefig('ligo-strain.png')

fft = ldata.fft()
print(fft)
plot = fft.abs().plot(xscale="log", yscale="log")
plot.show(warn=False)
plt.savefig('fft.png')

from scipy.signal import get_window
window = get_window('hann', ldata.size)
lwin = ldata * window
fftamp = lwin.fft().abs()
plot = fftamp.plot(xscale="log", yscale="log")
plot.show(warn=False)
plt.savefig('fft-window.png')

# asd = ldata.asd(fftlength=2, method="median")
# plot = asd.plot()
# plot.show(warn=False)

# ax = plot.gca()
# ax.set(xlim=(10, 1400), ylim=(1e-24, 1e-20))
# plot

# ldata2 = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512, cache=True)
# print("GW231123 data")
# print(ldata)
# print("")
lasd = ldata.asd(fftlength=4, method="median")
plot = lasd.plot()
ax = plot.gca()
ax.set_xlim(10, 1400)
ax.set_ylim(1e-24, 1e-20)
plot.show(warn=False)
# get Hanford data
hdata2 = TimeSeries.fetch_open_data('H1', int(gps)-512, int(gps)+512, cache=True)
hasd2 = hdata2.asd(fftlength=4, method="median")
# get Virgo data
# vdata2 = TimeSeries.fetch_open_data('V1', int(gps)-512, int(gps)+512, cache=True)
# vasd2 = vdata2.asd(fftlength=4, method="median")
# and plot using standard colours
ax.plot(hasd2, label='LIGO-Hanford', color='gwpy:ligo-hanford')
# ax.plot(vasd2, label='Virgo', color='gwpy:virgo')
# update the Livingston line to use standard colour, and have a label
lline = ax.lines[0]
lline.set_color('gwpy:ligo-livingston')  # change colour of Livingston data
lline.set_label('LIGO-Livingston')
ax.set_ylabel(r'Strain noise [$1/\sqrt{\mathrm{Hz}}$]')
ax.legend()
plot
plt.savefig('amplitude_spectral_density.png')

specgram = ldata.spectrogram2(fftlength=4, overlap=2, window='hann') ** (1/2.)
plot = specgram.plot()
plt.savefig('spec-no_filter.png')

ax = plot.gca()
ax.set_yscale('log')
ax.set_ylim(10, 1400)
ax.colorbar(
    clim=(1e-24, 1e-19),
    norm="log",
    label=r"Strain noise [$1/\sqrt{\mathrm{Hz}}$]",
)
plot
plt.savefig('spec-noise_zoom.png')

lq = ldata.q_transform(frange=(30, 500), qrange=(100, 110))
plot = lq.plot()
ax = plot.gca()
ax.set_epoch(gps)
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")
plt.savefig('spec-q_transform.png')

### Wave form

from pycbc.waveform import get_td_waveform, fd_approximants, td_approximants
# print('Time domain waveforms: ', td_approximants())
# print('Frequency domain waveforms: ', fd_approximants())

from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass

merger = Merger('GW231123_135430-v2')

strain = merger.strain('H1')

strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.savefig('waveform.png')
# plt.show()

conditioned = strain.crop(_,_)

plt.plot(conditioned.sample_times, conditioned)
plt.xlabel('Time(s)')

from pycbc.psd import interpolate, inverse_spectrum_truncation

psd = conditioned.psd(4)

psd = interpolate(psd, conditioned.delta_f)

psd = inverse_spectrum_truncation(psd, int(4 * conditioned_sample_rate),
                                  low_frequency_cutoff=15)







# ""Signal Model"""

m = 223
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                         mass1=m,
                         mass2=m,
                         delta_t=conditioned.delta_t,
                         f_lower=100)
hp.resize(len(conditioned))

# """to apporximate merger lcoation, change SNR time series to shift data into approximately first bin.

# cyclic_time_shift method.
# """

plt.figure()
plt.title('Before shifting')
plt.plot(hp.sample_times, hp)
plt.xlabel('Time (s)')
plt.ylabel('Strain')

template = hp.cyclic_time_shift(hp.start_time)

plt.figure()
plt.title('After shifting')
plt.plot(template.sample_times, template)
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.show()

# """Calculate Signal-to-noise time series"""

from pycbc.filter import matched_filter
import numpy

snr = matched_filter(template, conditioned,
                     psd=psd, low_frequency_cutoff=100)

# """upon generating waveform, remove corrupted time by the template filter and the psd filter (i.e. 4 sec at beginning or end). Longer signal would require more padding at beginning of vector."""

snr = snr.crop(4 +4, 4)

plt.figure(figsize=[10, 4])
plt.plot(snr.sample_times, abs(snr))
plt.ylabel('Signal-to-noise')
plt.xlabel('Time(s)')
plt.show

peak = abs(snr).numpy().argmax()
snrp = snr[peak]
time = snr.samplt_times[peak]

print("We found signal at {}s with SNR {}".format(time, abs(snrp)))

from pycbc.filter import sigma

dt = time - conditioned.start_time
aligned = template.cyclic_time_shift(dt)

aligned /= sigma(aligned, psd=psd, low_frequency_cutoff=100.0)

aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
aligned.start_time = conditioned.start_time

# mport matplotlib.pyplot as plt


m = Merger('GW231123_135430-v2')

ifos = ['H1', 'L1', 'V1']
data = {}
psd = {}

plt.figure(figsize=[20, 10])

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
plt.ylim(1e-47, 1e-41)
plt.xlim(20, 1024)
plt.ylabel('$Strain^2 / Hz$')
plt.xlabel('Frequency (Hz)')
plt.grid()
plt.legend()
plt.show()

from pycbc.waveform import get_fd_waveform
# from pycbc.filter import matched_filter

cmass = (m.median1d("mass1")+m.median1d("mass2")) / 2
cmass *= (1 + m.median1d("redshift"))

hp, _ =get_fd_waveform(approximant="IMRPhenomD",
                       mass1=cmass,
                       mass2=cmass,
                       f_lower=20.0,
                       delta_f=data[ifo].delta_f)
hp.resize(len(psd[ifo]))

snr={}
for ifo in ifos:
  snr[ifo] = matched_filter(hp, data[ifo], psd=psd[ifo], low_frequency_cutoff=100)
  snr[ifo] = snr[ifo].crop[10,9]

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