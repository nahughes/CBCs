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

########################################################################
########################################################################
########################################################################
########################################################################
############################ Wave form #################################
########################################################################
########################################################################
########################################################################
########################################################################


from pycbc.waveform import get_td_waveform, get_fd_waveform, fd_approximants, td_approximants
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


#See line 250 or 262 of Via_work to continue

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
plot.savefig('grid.png')