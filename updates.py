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

merger = Merger("GW231123")
strain = merger.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.show()

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

merger = Merger("GW231123")
strain = merger.strain('H1')
strain = highpass(strain, 15.0)
strain = resample_to_delta_t(strain, 1.0/2048)

plt.plot(strain.sample_times, strain)
plt.xlabel('Time(s)')
plt.show()

""" #Filter Wraparound"""
#may not be necessary to do with our target