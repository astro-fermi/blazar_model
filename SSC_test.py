################################################################################
# reordering and cleaning SSC_test for the upload on github
# testing the SSC_class
# 19 May 2017
# author: Cosimo Nigro (cosimonigro2@gmail.com)
################################################################################

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from ssc_model import model, numerics
import astropy.units as u
import timeit

time_grid = dict(time_min = 0., time_max = 3., time_bins = 50)
gamma_grid = dict(gamma_min = 2., gamma_max = 2e5, gamma_bins = 20)
emission_region = dict(R = 1e16, B = 1., t_esc = 1.5)
injected_spectrum = dict(norm = 1e-3, alpha = -2., t_inj = 1.5)
dist = 2*u.Mpc

start = timeit.default_timer()

# let us initialize the ssc object
SSC = model(time_grid, gamma_grid, emission_region, injected_spectrum)
num = numerics(SSC)

# let us evolve it
N_e = num.evolve()

# calculate the SED
# fetch the naima inverse compton and synchrotron object
SYN = num.synchrotron(N_e)
IC = num.inverse_compton(N_e)

stop = timeit.default_timer()
print 'Computational time: '
print stop - start, ' s'

# plotting section
fig, axes = plt.subplots(2, 1)
fig.subplots_adjust(hspace = 0.4)
font = {'family': 'serif',  'color':'black', 'weight': 'normal', 'size': 16.} # font definitions

# first plot with electron spectrum
axes[0].plot(SSC.gamma_grid, SSC.gamma_grid**2*SSC.N_e_inj,  ls = '--', lw=2, marker = '',
            color = 'turquoise', label = 'Injected Spectrum')
axes[0].plot(SSC.gamma_grid, SSC.gamma_grid**2*N_e,  ls = '-', lw=2, marker = '',
            color = 'crimson', label = 'numerical solution')
axes[0].legend(loc = 0, numpoints = 1., prop = {'size':12.})
axes[0].set_xscale('log')
axes[0].set_xlabel(r'$\gamma$')
axes[0].set_ylabel(r'$\gamma^2 \times n_{e}$')
axes[0].set_yscale('log')


energy = np.logspace(-7, 15, 100) * u.eV
axes[1].plot(energy, SYN.sed(energy, dist) + IC.sed(energy,dist), lw=3, color='royalblue', label='Synchrotron + Inverse Compton')
axes[1].set_xlabel(r'$E\,[eV]$')
axes[1].set_ylabel(r'$E^{2} \times {\rm d}F/{\rm d}E\,[erg\,cm^{-2}\,s^{-1}]$')
axes[1].set_ylim(1e-25, 1e-12)
axes[1].set_xscale('log')
axes[1].set_yscale('log')

fig.savefig('SSC_test_output.png')
plt.show()
