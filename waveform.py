########################################################################
########################################################################
########################################################################
########################################################################
############################ Wave form #################################
########################################################################
########################################################################
########################################################################
########################################################################

## Here I just worked through finding the data using pycbc instead of gwosc


# from pycbc.waveform import get_td_waveform, get_fd_waveform, fd_approximants, td_approximants
# print('Time domain waveforms: ', td_approximants())
# print('Frequency domain waveforms: ', fd_approximants())

# from pycbc.catalog import Merger
# from pycbc.filter import resample_to_delta_t, highpass 

import pycbc.catalog

# c = pycbc.catalog.Catalog(source='gwtc-4')
c= pycbc.catalog.find_event_in_catalog('GW231123_135430-v2')
print(c)
# print(c.names)