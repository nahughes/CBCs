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

# import matplotlib.pyplot as plt
# from pycbc.waveform.waveform import get_td_waveform, get_fd_waveform, fd_approximants, td_approximants
# print('Time domain waveforms: ', td_approximants())
# print("")
# print('Frequency domain waveforms: ', fd_approximants())
# print("")

# from pycbc.catalog import catalog, find_event_in_catalog, Merger
# from pycbc.filter import resample_to_delta_t, highpass 

import pycbc.catalog
event = pycbc.catalog.find_event_in_catalog('GW231123_135430-v2')
print(event)
print(f"Name: {event["commonName"]}")
merger = pycbc.catalog.Merger('GW231123_135430-v2')
mass1 = merger.strain('H1')
# print(strain)
# print(f"Name: {event["GPSstart"]}")

# listMergers = pycbc.catalog.Catalog(source='GWTC-4.0')
# # # print(c.names)
# substring_to_find = "231123"

# found_items = []
# for item in listMergers:
#     if substring_to_find in item:  # Checks if "an" is a substring of the item
#         found_items.append(item)

# if any(substring_to_find in item for item in listMergers):
#     print(f"At least one item contains '{substring_to_find}': {found_items}")
# else:
#     print(f"No item contains '{substring_to_find}'.")

# target = pycbc.catalog.Merger('GW231123_135430-v2')
# strain = target.strain('L1')
# print(strain)
# strain = highpass(strain, 15.0)
# strain = resample_to_delta_t(strain, 1.0/2048)

# plt.plot(strain.sample_times, strain)
# plt.xlabel('Time(s)')


'''
# The outputs of this function are the "plus" and "cross" polarizations
# of the gravitational-wave signal as viewed from the line of sight at
# a given source inclination (assumed face-on, i.e. zero inclination
# if not provided)
hp, hc = get_td_waveform(approximant="IMRPhenomD",
                         mass1=10,
                         mass2=10,
                         delta_t=1.0/16384,
                         f_lower=30)

plt.figure(figsize=plt.figaspect(0.4))
plt.plot(hp.sample_times, hp, label='Plus Polarization')
plt.plot(hp.sample_times, hc, label='Cross Polarization')
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.legend()
plt.grid()
plt.savefig('test1.png')

# Zoom in near the merger time
plt.figure(figsize=plt.figaspect(0.4))
plt.plot(hp.sample_times, hp, label='Plus Polarization')
plt.plot(hp.sample_times, hc, label='Cross Polarization')
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.xlim(-.01, .01)
plt.legend()
plt.grid()
plt.savefig('test2.png')
'''

'''



merger = Merger('GW231123')
strain = merger.strain('L1')
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

# import pycbc.catalog

# # c = pycbc.catalog.Catalog(source='gwtc-4')
# c= pycbc.catalog.find_event_in_catalog('GW231123_135430-v2')
# print(c)
# print(c.names)

########################################################################
########################################################################
########################################################################
########################################################################
############################ Wave form #################################
########################################################################
########################################################################
########################################################################
########################################################################


# from pycbc.waveform import get_td_waveform, get_fd_waveform, fd_approximants, td_approximants
# print('Time domain waveforms: ', td_approximants())
# print('Frequency domain waveforms: ', fd_approximants())

# from pycbc.catalog import  Merger
# c = pycbc.catalog.find_event_in_catalog('GW231123_135430')
# print(c)
# from pycbc.filter import resample_to_delta_t, highpass 

# import pycbc.catalg
import pycbc
merger = pycbc.catalog.Merger('GW231123_135430')

# import pycbc.catalog

# # c = pycbc.catalog.Catalog(source='gwtc-4')
# c= pycbc.catalog.find_event_in_catalog('GW231123_135430-v2')
# print(c)
# print(c.names)

'''