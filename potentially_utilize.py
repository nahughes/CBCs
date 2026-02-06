"""
Gravitational Wave Data Analysis
Event: GW231123_135430-v2
LIGO Detectors: L1 (Livingston), H1 (Hartford)
"""

import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

import numpy as np      #issues with numpy install :')
import matplotlib.pyplot as plt

# SECTION 1: DATA ACQUISITION AND EXPLORATION
def setup_gwosc():
    # pip install -q 'gwosc==0.7.1'
    import gwosc
    print(f"GWOSC Version: {gwosc.__version__}")
    return gwosc


def explore_datasets():
    from gwosc.datasets import find_datasets, event_gps
    from gwosc import datasets

    print("Current list of available catalogs: ")
    print(find_datasets(type="catalog"))

    #for GTWC-4 EVENTS
    gtwc4 = datasets.find_datasets(type='events', catalog='GTWC-4-confident')
    print(f'GWTC-4 events: {gtwc4}')

    #get GPS time for target event
    gps = event_gps('GW231123_135430-v2')
    print(f"GW231123_135430-v2 GPS: {gps}")

    #verify event at GPS time
    event = datasets.event_at_gps(1384782888)
    print(f"Event at GPS: {event}")

    #get run segment
    run_seg = datasets.run_segment('04')
    print(f"04 Run Segment: {run_seg}")

    return gps

#for general URLs
def get_event_data_urls(gps):
    from gwosc.locate import get_event_urls

    urls = get_event_urls('GW231123_135430-v2')
    print("Event URLs: ")
    print(urls)

    #Specify desired detector and duration
    urls_h1 = get_event_urls('GW231123_135430-v2', duration=32, detector='H1')
    print("\nH1 URLs (32sec duration): ")
    print(urls_h1)

    return urls

# SECTION 2: TIME SERIES ANALYSIS AND SPECTRAL ESTIMATION

def fetch_and_plot_timeseries(gps):
    # pip install -q 'gwpy==3.0.12'
    import gwpy
    from gwpy.timeseries import TimeSeries

    print(f"GWpy version: {gwpy.__version__}")

    #Define segment around selected event
    segment = (int(gps) - 10, int(gps) + 10)
    print(f"Analysis segment: {segment}")

    #Retrieve LIGO LIVINGSTON data
    ldata = TimeSeries.fetch_open_data('L1', *segment, verbose=True, cache=True)
    print(f"Livingston data: {ldata}")

    #Plot raw time series data
    plot = ldata.plot()
    plot()

    return ldata, segment

def compute_fft_analysis(ldata):
    from scipy.signal import get_window

    fft = ldata.fft()
    print(f"FFT: {fft}")

    #calculate PSD utilizing Hann Window
    window = get_window('hann', ldata.size)
    lwin = ldata * window
    fftamp = lwin.fft().abs()

    #plot fast fourier transform (FFT) amplitude
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.loglog(fft.frequencies.value, fft.abs().value)
    ax.set_xlabel('Frequency[Hz]')
    ax.set_ylabel('Amplitude')
    ax.set_title('Fast Fourier Transform Amplitude Spectrum')
    ax.grid(True, alpha=0.3)
    plt.savefig('FFT_plot.png', dpi=200) 
    plt.close()

    return fft

def compute_asd(ldata, gps):
    asd = ldata.asd(fftlength=2, method="median") 
    #compute amplitude spectral density (ASD) using Welch's method
    #reduce random fluctuations by averaging multiple FFT estimates (smoothes out)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.loglog(asd.frequencies.value, asd.value)
    ax.set_xlim(10, 1400)
    ax.set_ylim(1e-24, 1e-19)
    ax.set_xlabel('Frequency[Hz]')
    ax.set_ylabel(r'Strain Noise [$1/\sqrt{\mathrm{Hz}}$]')
    ax.set.title('Amplitude Spectral Density (ASD)')
    ax.grid(True, alpha=0.3)
    plt.savefig('ASD_plot.png', dpi=200) 
    plt.close()


    #optional: fetch longer segment for better averaging
    ldata2 = TimeSeries.fetch_open_data('L1', int(gps) - 512, int(gps) + 512, cache = True)
    lasd2 = ldata2.asd(fftlength=4, method="median")

    return lasd2, ldata2

def compare_detector_asds(gps):
    from gwpy.timeseries import TimeSeries #might not need double import

    ldata = TimeSeries.fetch_open_data('L1', int(gps) - 512, int(gps) + 512, cache = True)
    hdata = TimeSeries.fetch_open_data('H1', int(gps) - 512, int(gps) + 512, cache = True)

    #compute ASDs
    lasd = ldata.asd(fftlength=4, method="median")
    hasd = hdata.asd(fftlength=4, method="median")

    #comparison plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.loglog(lasd.frequencies.value, lasd.value,
              label='LIGO-Livingston', color='red')
    ax.loglog(hasd.frequencies.value, hasd.value,
              label='LIGO-Hanford', color='blue')
    ax.set_xlim(10, 1400)
    ax.set_ylim(1e-24, 1e-19)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'Strain Noise [$1/\sqrt{\mathrm{Hz}}$]')
    ax.set_title('Detector Noise Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('Detector_Noise_plot.png', dpi=200) 
    plt.close()

    return lasd, hasd



# SECTION 3: TIME-FREQUENCY ANALYSIS
def create_spectrogram(gps):
    from gwpy.timeseries import TimeSeries #might not need import, here in case

    ldata = TimeSeries.fetch_open_data('L1', int(gps) - 512, int(gps) + 512, cache=True)

    specgram = ldata.spectrogram2(fftlength=4, overlap=2, window='hann') ** 0.5

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh(specgram.times.value, specgram.frequencies.value,
                       specgram.value, cmap='coolwarm', shading='auto')     
    ax.set_yscale('log')
    ax.set_ylim(10, 1400)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Frequency [Hz]')
    ax.set_title('Spectrogram')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'Strain noise [$1/\sqrt{\mathrm{Hz}}$]')
    im.set_clim(1e-24, 1e-19)

    plt.savefig('Time_series_spectrogram.png', dpi=200)
    plt.close()

    return specgram

def q_transform_analysis(gps):
    from gwpy.timeseries import TimeSeries #might not need import, here in case
    
    #fetch data segment around merger
    segment = (int(gps) - 30, int(gps) + 20)
    hdata = TimeSeries.fetch_open_data('H1', *segment, verbose=True, cache=True)
    
    #Compute Q-transform
    hq = hdata.q_transform(frange=(30, 500))
    
    #Plot Q-transform
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh(hq.times.value, hq.frequencies.value, 
                       hq.value, cmap='viridis', shading='auto')
    ax.set_epoch(gps)
    ax.set_yscale('log')
    ax.set_ylim(30, 500)
    ax.set_xlabel('Time from merger [s]')
    ax.set_ylabel('Frequency [Hz]')
    ax.set_title('Q-transform')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Normalized energy')

    plt.savefig('Q_transform_plot.png', dpi=200)
    plt.close()

    #zoomed in Q-transform around merger
    hq_zoom = hdata.q_transform(frange=(25, 100), qrange=(80, 110), 
                                 outseg=(gps - 0.8, gps + 0.8))
    
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh(hq_zoom.times.value, hq_zoom.frequencies.value, 
                       hq_zoom.value, cmap='coolwarm', shading='auto')
    ax.set_epoch(gps)
    ax.set_yscale('log')
    ax.set_xlabel('Time from merger [s]')
    ax.set_ylabel('Frequency [Hz]')
    ax.set_title('Q-transform (Zoomed around merger)')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Normalized energy')
    
    plt.savefig('Zoomed_In_Q_transform_plot.png', dpi=200)
    plt.close()
    
    return hq

def apply_gating(gps):
    from gwpy.timeseries import TimeSeries

    segment = (int(gps) - 30, int(gps) + 20)    
    ldata = TimeSeries.fetch_open_data('L1', *segment, verbose=True, cache=True)

    #Apply gate to remove glitch
    gated_ldata = ldata.gate(tzero=0.25, tpad=0.25)
    
    #Compute Q-transform of gated data
    gated_lq = gated_ldata.q_transform(frange=(30, 500), qrange=(100, 110))
    
    #Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh(gated_lq.times.value, gated_lq.frequencies.value, 
                       gated_lq.value, cmap='viridis', shading='auto')
    ax.set_epoch(gps)
    ax.set_yscale('log')
    ax.set_xlabel('Time from merger [s]')
    ax.set_ylabel('Frequency [Hz]')
    ax.set_title('Q-transform (After Gating)')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Normalized energy')
    im.set_clim(0, 20)
    
    plt.savefig('Gated_Ldata_plot.png', dpi=200)
    plt.close()
    
    return gated_ldata


# SECTION 4: WAVEFORM GENERATION (APPLIED AND THEORETICAL)
    #Generate theoretical gravitational waveforms for different masses
    #Utilize PyCBC to generate waveforms with IMRPhenomD approximant

def generate_waveforms(sample_rate=4096):
    # pip install -q PyCBC==2.4.1 lalsuite==7.25
    from pycbc.waveform import get_td_waveform, td_approximants, fd_approximants
    
    print('Time domain waveforms:', td_approximants())
    print('Frequency domain waveforms:', fd_approximants())
    
    #Generate waveforms for different total masses
    #GW231123: m1=137, m2=101 M_sun (total ~238 M_sun)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for m in [101, 137, 223]:
        hp, hc = get_td_waveform(
            approximant="IMRPhenomD",
            mass1=m,
            mass2=m,
            delta_t=1.0/sample_rate,
            f_lower=20  # Low frequency cutoff
        )
        ax.plot(hp.sample_times, hp, label=f'$M_{{1,2}}={m}$ $M_\odot$')
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Strain')
    ax.set_title('Gravitational Waveforms for Different Masses')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('Waveform_231123_TD_FD.png', dpi=200)
    plt.close()
    
    return hp

# SECTION 5: MATCHED FILTERING

#Demonstrate matched filtering for GW231123 (optimal detection technique in Gaussian noise)

def matched_filtering_example(gps, sample_rate=4096):
    # pip install -q 'lalsuite==7.25' 'PyCBC==2.6.0'
    from pycbc.waveform import get_td_waveform
    from pycbc.filter import matched_filter
    from pycbc.catalog import Merger
    from pycbc.filter import resample_to_delta_t, highpass
    from pycbc.psd import interpolate, inverse_spectrum_truncation
    
    #Load desired event data
    merger = Merger('GW231123_135430-v2')
    strain = merger.strain('H1')
    
    #Precondition: highpass filter and resample
    strain = highpass(strain, 15.0)
    strain = resample_to_delta_t(strain, 1.0/2048)

    #plot preconditioned strain
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(strain.sample_times, strain)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Strain')
    ax.set_title('Preconditioned Strain Data (H1)')
    plt.savefig('Preconditioned_strain.png', dpi=200)
    plt.close()

    #remove filter wraparound
    conditioned = strain.crop(2, 2)

    #estimate the PSD using Welch's method
    psd = conditioned.psd(4)
    psd = interpolate(psd, conditioned.delta_f)
    psd = inverse_spectrum_truncation(
        psd, 
        int(4 * conditioned.sample_rate),
        low_frequency_cutoff=15
    )

    #generate template waveform (use component masses for target)
    #using GW231123 in this instance

    m1, m2 = 137, 101 #solar masses
    hp, hc = get_td_waveform(
        approximant="SEOBNRv4_opt",
        mass1=m1,
        mass2=m2,
        delta_t=conditioned.delta_t,
        f_lower=20
    )

    #resize the template to match data length
    hp.resize(len(conditioned))
    
    #Shift template to start at t=0
    template = hp.cyclic_time_shift(hp.start_time)
    
    #Compute SNR time series
    snr = matched_filter(template, conditioned, psd=psd, low_frequency_cutoff=20)
    
    #Crop extra corrupted data at the edges
    snr = snr.crop(4 + 4, 4)
    
    #Plot SNR time series
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(snr.sample_times, abs(snr))
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Signal-to-Noise Ratio')
    ax.set_title('Matched Filter SNR Time Series')
    ax.grid(True, alpha=0.3)
    plt.savefig('SNR_timeseries.png', dpi=200)
    plt.close()
    
    #Find peak SNR
    peak = abs(snr).numpy().argmax()
    snrp = snr[peak]
    time = snr.sample_times[peak]
    
    print(f"Found signal at {time}s with SNR {abs(snrp):.2f}")
    
    return snr, template, conditioned, psd

def visualize_matched_signal(merger_time, template, conditioned, psd, snr):
    from pycbc.filter import sigma

    #get peak SNR value
    peak = abs(snr).numpy().argmax()
    snrp = snr[peak]
    time = snr.sample_times[peak]
    
    #Align the template to detected time
    dt = time - conditioned.start_time
    aligned = template.cyclic_time_shift(dt)
    
    #Normalize template
    aligned /= sigma(aligned, psd=psd, low_frequency_cutoff=20.0)
    aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
    aligned.start_time = conditioned.start_time
    
    #Whiten both the data and template
    white_data = (conditioned.to_frequencyseries() / psd**0.5).to_timeseries()
    white_template = (aligned.to_frequencyseries() / psd**0.5).to_timeseries()
    
    #Bandpass btwn 30-300 Hz
    white_data = white_data.highpass_fir(30., 512).lowpass_fir(300, 512)
    white_template = white_template.highpass_fir(30, 512).lowpass_fir(300, 512)
    
    #Select time slice around merger
    white_data = white_data.time_slice(merger_time - 0.2, merger_time + 0.1)
    white_template = white_template.time_slice(merger_time - 0.2, merger_time + 0.1)
    
    #Plot overlay
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.plot(white_data.sample_times, white_data, label="Data", alpha=0.7)
    ax.plot(white_template.sample_times, white_template, label="Template", alpha=0.7)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Whitened Strain')
    ax.set_title('Whitened Data and Template Overlay')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('whitened_overlay.png', dpi=200)
    plt.close()
    
    #Subtract the signal from data
    subtracted = conditioned - aligned
    
    #Create Q-transforms before and after subtraction
    fig, axes = plt.subplots(2, 1, figsize=(15, 8))
    
    for idx, (data, title) in enumerate([
        (conditioned, 'Original H1 Data'),
        (subtracted, 'Signal Subtracted from H1 Data')
    ]):
        t, f, p = data.whiten(4, 4).qtransform(
            0.001, logfsteps=100, qrange=(8, 8), frange=(20, 512)
        )
        im = axes[idx].pcolormesh(t, f, p**0.5, vmin=1, vmax=6, 
                                  cmap='viridis', shading='auto')
        axes[idx].set_yscale('log')
        axes[idx].set_xlabel('Time [s]')
        axes[idx].set_ylabel('Frequency [Hz]')
        axes[idx].set_title(title)
        axes[idx].set_xlim(merger_time - 2, merger_time + 1)
        plt.colorbar(im, ax=axes[idx])
    
    plt.tight_layout()
    plt.savefig('Signal_subtraction.png', dpi=200)
    plt.close()
    
    return aligned

# SECTION 6: SIGNAL CONSISTENCY AND SIGNIFICANCE
def multi_detector_analysis(gps):
    # pip install -q 'PyCBC==2.6.0'
    from pycbc.catalog import Merger
    from pycbc.filter import resample_to_delta_t, highpass, matched_filter
    from pycbc.psd import interpolate, inverse_spectrum_truncation
    from pycbc.waveform import get_fd_waveform
    from pycbc.vetoes import power_chisq
    from pycbc.events.ranking import newsnr
    from pycbc.detector import Detector
    
    m = Merger("GW231123_135430-v2")
    ifos = ['H1', 'L1']
    data = {}
    psd = {}
    
    #Process each detector
    for ifo in ifos:
        ts = m.strain(ifo).highpass_fir(15, 512)
        data[ifo] = resample_to_delta_t(ts, 1.0/2048).crop(2, 2)
        
        p = data[ifo].psd(2)
        p = interpolate(p, data[ifo].delta_f)
        p = inverse_spectrum_truncation(
            p, 
            int(2 * data[ifo].sample_rate), 
            low_frequency_cutoff=15.0
        )
        psd[ifo] = p
    
    #Plot the PSDs
    fig, ax = plt.subplots(figsize=(10, 6))
    for ifo in ifos:
        ax.loglog(psd[ifo].sample_frequencies, psd[ifo], label=ifo)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-47, 1e-43)
    ax.set_xlim(20, 1024)
    ax.set_ylabel('$Strain^2 / Hz$')
    ax.set_xlabel('Frequency [Hz]')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('detector_PSDs.png', dpi=200)
    plt.close()
    
    #Generate a waveform
    cmass = (m.median1d("mass1") + m.median1d("mass2")) / 2
    cmass *= (1 + m.median1d("redshift"))
    
    hp, _ = get_fd_waveform(
        approximant="IMRPhenomD",
        mass1=cmass, 
        mass2=cmass,
        f_lower=20.0, 
        delta_f=data['H1'].delta_f
    )
    hp.resize(len(psd['H1']))
    
    #Compute the SNR for each detector
    snr = {}
    for ifo in ifos:
        snr[ifo] = matched_filter(hp, data[ifo], psd=psd[ifo], low_frequency_cutoff=20)
        snr[ifo] = snr[ifo].crop(5, 4)
    
    #Plot SNR time series
    fig, ax = plt.subplots(figsize=(14, 4))
    for ifo in ifos:
        ax.plot(snr[ifo].sample_times, abs(snr[ifo]), label=ifo)
    ax.set_xlim(m.time - 0.15, m.time + 0.15)
    ax.set_ylim(0, 15)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Signal-to-Noise Ratio')
    ax.set_title('Multi-Detector SNR (Zoomed)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('multi_detector_SNR.png', dpi=200)
    plt.close()
    
    #Compute chi-squared statistic
    chisq = {}
    nbins = 30
    dof = nbins * 2 - 2
    
    for ifo in ifos:
        chisq[ifo] = power_chisq(hp, data[ifo], nbins, psd[ifo], 
                                 low_frequency_cutoff=20.0)
        chisq[ifo] = chisq[ifo].crop(5, 4)
        chisq[ifo] /= dof
    
    #Plot chi-squared
    fig, ax = plt.subplots(figsize=(14, 4))
    for ifo in ifos:
        ax.plot(chisq[ifo].sample_times, chisq[ifo], label=ifo)
    ax.set_xlim(m.time - 0.15, m.time + 0.15)
    ax.set_ylim(0, 5)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('$\\chi^2_r$')
    ax.set_title('Chi-squared Test Statistic')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('Chi_squared.png', dpi=200)
    plt.close()
    
    #Compute re-weighted SNR
    nsnr = {ifo: newsnr(abs(snr[ifo]), chisq[ifo]) for ifo in ifos}
    
    #Plot re-weighted SNR
    fig, ax = plt.subplots(figsize=(14, 4))
    for ifo in ifos:
        ax.plot(snr[ifo].sample_times, nsnr[ifo], label=ifo)
    ax.set_xlim(m.time - 0.15, m.time + 0.15)
    ax.set_ylim(0, 15)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Re-weighted SNR')
    ax.set_title('Re-weighted SNR Time Series')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('reweighted_SNR.png', dpi=200)
    plt.close()
    
    return snr, chisq, nsnr


def calculate_significance(m, nsnr, snr):
    from pycbc.detector import Detector
    
    ifos = ['H1', 'L1']
    
    #Calculate light-travel time between detectors
    d = Detector("H1")
    tof = {'L1': d.light_travel_time_to_detector(Detector("L1"))}
    
    #Find peak times
    ptime = {}
    for ifo in ifos:
        if ifo != 'H1':
            ptime[ifo] = snr[ifo].sample_times[nsnr[ifo].argmax()]
    
    #Define coincidence window
    #Note: Need to calculate proper window for signal timing
    print("\nSignificance Analysis:")
    print(f"H1 peak SNR: {nsnr['H1'].max():.2f}")
    print(f"L1 peak SNR: {nsnr['L1'].max():.2f}")
    print(f"Light travel time H1-L1: {tof['L1']*1000:.2f} ms")
    
    return tof

# SECTION 7: PARAMETER ESTIMATION W/ BILBY
def setup_bilby_analysis(gps, time_of_event):
    
    #Set up Bilby parameter estimation analysis.
    #Bilby performs Bayesian inference to estimate source parameters.

    #pip install -U -q bilby==2.6.0 dynesty==2.1.5 corner==2.2.3
    import bilby
    from bilby.core.prior import Uniform, PowerLaw
    from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters
    from gwpy.timeseries import TimeSeries
    
    bilby.core.utils.log.setup_logger(log_level="WARNING")
    
    print(f"Bilby version: {bilby.__version__}")
    
    #Initialize interferometers
    H1 = bilby.gw.detector.get_empty_interferometer("H1")
    L1 = bilby.gw.detector.get_empty_interferometer("L1")
    
    #Define analysis parameters
    post_trigger_duration = 2
    duration = 4
    analysis_start = time_of_event + post_trigger_duration - duration
    
    #Fetch strain data
    H1_data = TimeSeries.fetch_open_data(
        "H1", analysis_start, analysis_start + duration, 
        sample_rate=4096, cache=True
    )
    L1_data = TimeSeries.fetch_open_data(
        "L1", analysis_start, analysis_start + duration, 
        sample_rate=4096, cache=True
    )
    
    #Set strain data
    H1.set_strain_data_from_gwpy_timeseries(H1_data)
    L1.set_strain_data_from_gwpy_timeseries(L1_data)
    
    #Fetch PSD data (longer duration for better estimate)
    psd_duration = duration * 32
    psd_start_time = analysis_start - psd_duration
    
    H1_psd_data = TimeSeries.fetch_open_data(
        "H1", psd_start_time, psd_start_time + psd_duration,
        sample_rate=4096, cache=True
    )
    L1_psd_data = TimeSeries.fetch_open_data(
        "L1", psd_start_time, psd_start_time + psd_duration,
        sample_rate=4096, cache=True
    )
    
    #Compute PSDs
    psd_alpha = 2 * H1.strain_data.roll_off / duration
    H1_psd = H1_psd_data.psd(
        fftlength=duration, overlap=0, 
        window=("tukey", psd_alpha), method="median"
    )
    L1_psd = L1_psd_data.psd(
        fftlength=duration, overlap=0, 
        window=("tukey", psd_alpha), method="median"
    )
    
    #Initialize PSDs
    H1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=H1_psd.frequencies.value, 
        psd_array=H1_psd.value
    )
    L1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=L1_psd.frequencies.value, 
        psd_array=L1_psd.value
    )
    
    #Set maximum frequency
    H1.maximum_frequency = 1024
    L1.maximum_frequency = 1024
    
    #Plot strain data and PSD
    fig, ax = plt.subplots(figsize=(10, 6))
    idxs = H1.strain_data.frequency_mask
    ax.loglog(H1.strain_data.frequency_array[idxs],
              np.abs(H1.strain_data.frequency_domain_strain[idxs]),
              label='Strain data', alpha=0.7)
    ax.loglog(H1.power_spectral_density.frequency_array[idxs],
              H1.power_spectral_density.asd_array[idxs],
              label='ASD', linewidth=2)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Strain [strain/$\\sqrt{Hz}$]')
    ax.set_title('H1 Strain Data and ASD')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('Bilby_Data_PSD.png', dpi=200)
    plt.close()
    
    return H1, L1, [H1, L1]

def define_priors(time_of_event):
    #Define prior distributions for parameter estimation
    import bilby #duplicate import, might need
    from bilby.core.prior import Uniform, PowerLaw
    
    #Create prior dictionary
    #Based on GW231123: m1~137, m2~101 M_sun
    prior = bilby.core.prior.PriorDict()
    
    #Mass parameters
    prior['chirp_mass'] = Uniform(
        name='chirp_mass', minimum=80.0, maximum=120.0
    )
    prior['mass_ratio'] = Uniform(
        name='mass_ratio', minimum=0.5, maximum=1.0
    )
    
    #Spin parameters (currently set to zero for simplicity)
    prior['a_1'] = 0.0
    prior['a_2'] = 0.0
    prior['tilt_1'] = 0.0
    prior['tilt_2'] = 0.0
    prior['phi_12'] = 0.0
    prior['phi_jl'] = 0.0
    
    #Extrinsic parameters
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(
        name="geocent_time", 
        minimum=time_of_event - 0.1, 
        maximum=time_of_event + 0.1
    )
    
    #Sky location (set to approximate values, should be fit)
    prior['dec'] = 0.0  #Placeholder: update with actual sky location
    prior['ra'] = 0.0   #Placeholder: update with actual sky location
    prior['theta_jn'] = Uniform(name='theta_jn', minimum=0, maximum=np.pi)
    prior['psi'] = Uniform(name='psi', minimum=0, maximum=np.pi)
    
    #Distance
    prior['luminosity_distance'] = PowerLaw(
        alpha=2, name='luminosity_distance', 
        minimum=50, maximum=5000, unit='Mpc'
    )
    
    return prior

def run_parameter_estimation(interferometers, prior, label="GW231123"):
    """
    Run parameter estimation using Bilby.
    This is computationally intensive and may take hours to complete.
    """
    import bilby
    from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters
    
    #Define waveform generator
    waveform_arguments = dict(
        waveform_approximant='IMRPhenomXP',
        reference_frequency=100.,
        catch_waveform_errors=True
    )
    
    waveform_generator = bilby.gw.WaveformGenerator(
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        waveform_arguments=waveform_arguments,
        parameter_conversion=convert_to_lal_binary_black_hole_parameters
    )
    
    #Create likelihood
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers, 
        waveform_generator, 
        priors=prior,
        time_marginalization=True,
        phase_marginalization=True,
        distance_marginalization=True
    )
    
    #Run sampler (simplified settings for demonstration)
    print("\nRunning parameter estimation...")
    print("Note: This is a demonstration with fast settings.")
    print("For production analysis, increase nlive and remove dlogz constraint.")
    
    result = bilby.run_sampler(
        likelihood, 
        prior, 
        sampler='dynesty',
        outdir='/mnt/user-data/outputs/pe_results',
        label=label,
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
        nlive=250,      #Low for speed - use 1000+ for real analysis
        dlogz=1.0,      #Fast termination - remove for real analysis
        clean=True,
    )
    
    return result


def analyze_pe_results(result):
    #Analyze and plot parameter estimation results.
    import corner
    
    #Print summary statistics
    print("\n=== Parameter Estimation Results ===")
    
    #Chirp mass
    Mc = result.posterior["chirp_mass"].values
    lower_bound = np.quantile(Mc, 0.05)
    upper_bound = np.quantile(Mc, 0.95)
    median = np.quantile(Mc, 0.5)
    
    print(f"\nChirp Mass:")
    print(f"  Median: {median:.2f} M_sun")
    print(f"  90% C.I.: [{lower_bound:.2f}, {upper_bound:.2f}] M_sun")
    
    #Plot chirp mass histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(result.posterior["chirp_mass"], bins=30, alpha=0.7, edgecolor='black')
    ax.axvspan(lower_bound, upper_bound, color='C1', alpha=0.3, 
               label='90% Credible Interval')
    ax.axvline(median, color='C1', linewidth=2, label='Median')
    ax.set_xlabel('Chirp Mass [$M_\\odot$]')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Chirp Mass Posterior Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('Chirp_Mass_Posterior.png', dpi=200)
    plt.close()
    
    #Create a corner plot
    result.plot_corner(
        parameters=["chirp_mass", "mass_ratio", "geocent_time", "luminosity_distance"],
        prior=True,
        save=True,
        filename='Corner_plot.png'
    )
    
    #Print the Bayes factor
    print(f"\nln(Bayes Factor): {result.log_bayes_factor:.2f} ± {result.log_evidence_err:.2f}")
    
    return result

# SECTION 8: WORKING W/ PUBLISHED POSTERIORS

def load_published_posteriors(event_name='GW231123'):
    #Load and analyze published posterior samples from GWTC.
    #Note!!!: Requires downloading the posterior file from GWTC.

    import h5py
    import pandas as pd
    import corner
    
    #File should be downloaded from GWTC-4 release
    posterior_file = f'./{event_name}_GWTC-4.hdf5'
    
    try:
        posterior = h5py.File(posterior_file, 'r')
        
        print("File datasets:", list(posterior.keys()))
        
        #Load overall posterior
        samples = pd.DataFrame.from_records(
            np.array(posterior['Overall_posterior'])
        )
        
        print("\nPosterior parameters:")
        print(samples.columns.tolist())
        
        return posterior, samples
        
    except FileNotFoundError:
        print(f"\nPosterior file not found: {posterior_file}")
        print("Please download from GWTC-4 data release:")
        print("https://zenodo.org/records/8177023") #DOUBLE CHECK THIS LINK
        return None, None


def plot_published_posteriors(posterior, samples):
    #Create visualizations of published posterior samples.
    import corner
    
    if posterior is None:
        return
    
    #Create a corner plot of key parameters
    params_to_plot = ['luminosity_distance_Mpc', 'm1_detector_frame_Msun', 
                      'm2_detector_frame_Msun']
    
    #Check which parameters are available
    available_params = [p for p in params_to_plot if p in samples.columns]
    
    if available_params:
        fig = corner.corner(
            samples[available_params].values,
            labels=[p.replace('_', ' ') for p in available_params],
            show_titles=True,
            title_fmt='.2f'
        )
        plt.savefig('Published_corner.png', dpi=200)
        plt.close()
    
    #Plot luminosity distance comparison
    if 'luminosity_distance_Mpc' in samples.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for key in posterior.keys():
            if 'posterior' in key.lower():
                try:
                    dist = posterior[key]['luminosity_distance_Mpc']
                    ax.hist(dist, bins=50, alpha=0.5, density=True, label=key)
                except:
                    pass
        
        ax.set_xlabel('Luminosity Distance [Mpc]')
        ax.set_ylabel('Probability Density')
        ax.set_title('Luminosity Distance Posteriors')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.savefig('Distance_comparison.png', dpi=200)
        plt.close()


def compute_source_frame_masses(samples):
    #Convert detector frame masses to source frame.
    #Accounts for cosmological redshift.

    from astropy.cosmology import Planck15, z_at_value
    import astropy.units as u
    
    if samples is None:
        return None
    
    #Compute redshift for each sample
    print("\nComputing source frame masses... ")
    
    z_values = []
    for dist in samples['luminosity_distance_Mpc'][:100]:  #Limit for speed
        z = z_at_value(Planck15.luminosity_distance, dist * u.Mpc)
        z_values.append(z.value)
    
    z_array = np.array(z_values)
    
    #Compute source frame masses
    m1_source = samples['m1_detector_frame_Msun'][:100] / (1 + z_array)
    m2_source = samples['m2_detector_frame_Msun'][:100] / (1 + z_array)
    
    #Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].hist2d(samples['m1_detector_frame_Msun'][:100],
                   samples['m2_detector_frame_Msun'][:100],
                   bins=20, cmap='Blues')
    axes[0].set_xlabel('$m_1$ (detector frame) [$M_\\odot$]')
    axes[0].set_ylabel('$m_2$ (detector frame) [$M_\\odot$]')
    axes[0].set_title('Detector Frame Masses')
    
    axes[1].hist2d(m1_source, m2_source, bins=20, cmap='Blues')
    axes[1].set_xlabel('$m_1$ (source frame) [$M_\\odot$]')
    axes[1].set_ylabel('$m_2$ (source frame) [$M_\\odot$]')
    axes[1].set_title('Source Frame Masses')
    
    plt.tight_layout()
    plt.savefig('Mass_frame_comparsion.png', dpi=200)
    plt.close()
    
    return z_array, m1_source, m2_source

# MAIN EXECUTION (TRIAL)

def main():
    #Main execution function.
    #Runs a complete analysis pipeline for GW231123.
    
    print("=" * 70)
    print("GW ANALYSIS: GW231123_135430-v2")
    print("=" * 70)
    
    #GPS time of event
    gps = 1384782888.6
    time_of_event = gps
    
    try:
        #Section 1: Data Acquisition
        print("\n--- SECTION 1: Data Acquisition ---")
        gwosc = setup_gwosc()
        gps = explore_datasets()
        urls = get_event_data_urls(gps)
        
        #Section 2: Time Series Analysis
        print("\n--- SECTION 2: Time Series Analysis ---")
        ldata, segment = fetch_and_plot_timeseries(gps)
        fft = compute_fft_analysis(ldata)
        lasd, ldata2 = compute_asd(ldata, gps)
        lasd_comp, hasd_comp = compare_detector_asds(gps)
        
        #Section 3: Time-Frequency Analysis
        print("\n--- SECTION 3: Time-Frequency Analysis ---")
        specgram = create_spectrogram(gps)
        hq = q_transform_analysis(gps)
        gated_data = apply_gating(gps)
        
        #Section 4: Waveform Generation
        print("\n--- SECTION 4: Waveform Generation ---")
        hp = generate_waveforms()
        
        #Section 5: Matched Filtering
        print("\n--- SECTION 5: Matched Filtering ---")
        snr, template, conditioned, psd = matched_filtering_example(gps)
        aligned = visualize_matched_signal(gps, template, conditioned, psd, snr)
        
        #Section 6: Multi-Detector Analysis
        print("\n--- SECTION 6: Multi-Detector Analysis ---")
        snr_multi, chisq, nsnr = multi_detector_analysis(gps)
        
        print("\n" + "=" * 70)
        print("Analysis complete! Output files saved to /mnt/user-data/outputs/")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()