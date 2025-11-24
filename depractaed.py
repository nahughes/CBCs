# Old version of Waveform from Friday Nov. 21


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