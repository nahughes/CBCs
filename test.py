import gwpopulation
import numpy as np
import bilby
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import os 

""" #Reading in the 03b Data """
#find HDF5 files in the Parameter Estimation Release: https://zenodo.org/records/17014085 <- record somewhere on here
# !mkdir tuto_A.2
# !wget --no-verbose -O ./tuto_A.2/downsampled.tar.gz "https://dcc.ligo.org/public/0182/T2200137/001/O3bPE_downsampled.tar.gz"
# !tar -xvzf ./tuto_A.2/downsampled.tar.gz -C ./tuto_A.2
# import glob
# files = glob.glob('./tuto_A.2/*.h5')

posteriors = []
for file in files:
    if '200105' in file or '200115' in file:
        continue
    with h5py.File(file, 'r') as eventfile:
        post = pd.DataFrame()
        samples = eventfile['C01:Mixed']['posterior_samples']
        post['mass_1'] = samples['mass_1_source']
        posteriors.append(post)

""" #Defining the Population Model """
#go over text w/ Dr. Hughes
population = gwpopulation.models.mass.power_law_primary_mass_ratio

""" #Selection Effects Term """
#go over text w/ Dr. Hughes
m1 = np.linspace(3.1, 60, 100)
q = np.linspace(0.05, 1, 100)
m1_grid, q_grid = np.meshgrid(m1, q)
chirps = bilby.gw.conversion.component_masses_to_chirp_mass(m1_grid, m1_grid * q_grid)
vt_selection = chirps**(5/2)

mass_grid = dict()
mass_grid['mass_1'] = m1_grid
mass_grid['mass_ratio'] = q_grid

def p_det(params):
    p_m = gw.population.models.mass.power_law_primary_mass_ratio(mass_grid, **params)

    return np.trapz(np.trapz(p_m * vt_selection, q, axis = 0), m1, axis=0)

""" #Setting up the Run """
likelihood = gwpopulation.hyperpe.HyperparameterLikelihood(posterios = posterios, hyper_prior = population, selection_function = p_det)

#define priors for the population hyperparameters (assume uniform priors)
from bilby.core.prior import PriorDict, Uniform

priors = PriorDict(
    dict(
        alpha = Uniform(minimum=0, maximum=4, latex_label='$\\alpha$'),
        beta = Uniform(minimum=0, maximum=7, name='beta', latex_label='$\\beta_{q}$'),
        mmax = Uniform(minimum=30, maximum=60, name='mmax', latex_label='$m_{\\max}$'),
        mmin = Uniform(minimum=3, maximum=12, name='mmin', latex_label='$m_{\\min}$')
))

""" #Running the Sampler """
bilby.run_sampler(likelihood = likelihood, priors = priors, label = 'O3b', use_ratio=True, resume=False, outdir = './tuto_A.2/O3b_result', sampler='nestle')

""" #Viewing our Results """
#use bilby to read in the result file and make a corner plot to see inferred posterior distribution
result = bilby.result.read_in_result('./tuto_A.2/O3b_result.json')
result.plot_corner(save=False)

#sample model to show how we may prepare GWPopulation analyses and view the results
#NOTE: the plot_corner() command both the corner plot inline and can save a copt in the directory we save our result file
#we now look at our inferred distribution for m1. Plot the inferred 90% C.I.

m1s = np.linspace(2, 70, 1000)
qs = np.linspace(0, 1, 1000)

mlines = []
samples = result.posterior.sample(1000)
for ii in range (len(samples)):
    sample = samples.iloc[ii]
    m1probs = gwpopulations.utils.powerlaw(m1s, alpha = -sample['alpha'], high=sample['mmax'], low=sample['mmin'])
    mlines.append(m1probs)
mlines = np.array(mlines)

plt.fill_between(m1s, np.percentile(mlines, 5, axis=0), np.percentile(mlines, 95, axis=0), alpha = 0.4)
plt.xlim(4, 60)
plt.yscale('log')
plt.ylabel(r'$p(m_1)$')
plt.xlabel(r'$m_1$')