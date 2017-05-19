from __future__ import division
from math import pi
import numpy as np
from constants import *

class model:
    '''
    This class creates an object which contains the particle spectra of electrons
    and photons in case of SSC process in AGN.
    Reference is Chiaberge & Ghisellini (1998) MNRAS, 306,551 (1999)
    '''

    def __init__(self, time_grid, gamma_grid, emission_region, injected_spectrum):
        '''
        This is the constructor of our SSC process with
        * time_grid = dict(minimum, maximum time, binning)
        initializing the time grid
        * gamma_grid = dict(minimum, maximum Lorentz factor, binning)
        initializing the Lorentz factor grid
        * emission_region dict(radius, magnetic field, escaping time)
        with parameters describing the emission region, escaping time in order of R/c
        * injected_spectrum dict(norm, index)
        initializing the injected spectrum (a powerlaw for now).
        '''

        # emission regions attributes definitions
        self.R = emission_region['R'] # in cm
        self.crossing_time = self.R/c
        self.B = emission_region['B'] # in Gauss
        self.U_B = self.B/(8.*pi) # magnetic field density
        self.t_esc = emission_region['t_esc']*self.crossing_time

        # time grid attributes definition
        self.time_min = time_grid['time_min']*self.crossing_time
        self.time_max =  time_grid['time_max']*self.crossing_time
        self.time_bins = time_grid['time_bins']
        self.delta_t = (self.time_max - self.time_min)/self.time_bins
        self.time_grid = np.linspace(self.time_min, self.time_max, self.time_bins)

        # gamma grid attributes defintion
        self.gamma_min = gamma_grid['gamma_min']
        self.gamma_max = gamma_grid['gamma_max']
        self.gamma_bins = gamma_grid['gamma_bins']
        # following Eq.(5) of the Reference, LorentzFactor grid has to be logspaced
        # we loop over range(-1,N+2) include \gamma_{-1} and \gamma{N+1}
        gamma_grid = []
        for j in range(-1, self.gamma_bins + 2):
            gamma_grid.append(self.gamma_min*(self.gamma_max/self.gamma_min)**((j-1)/(self.gamma_bins-1)))
        # gamma_grid_ext will be the grid that contains \gamma_{-1} and \gamma{N+1}
        self.gamma_grid_ext = np.array(gamma_grid)
        # gamma_grid will be the grid with elements from \gamma_{0} to \gamma_{N}
        self.gamma_grid = self.gamma_grid_ext[1:-1]
        # energy grid for use in calculating the synchrotron emission
        self.energy_grid = self.gamma_grid*E_rest

        # let's create an array for the density of radiation
        self.U_rad = np.zeros(len(self.gamma_grid), float)

        # injected spectrum attributes
        self.inj_spectr_norm = injected_spectrum['norm']
        self.inj_spectr_index = injected_spectrum['alpha']
        self.inj_time = injected_spectrum['t_inj']*self.crossing_time # injection time


    @property
    def gamma_grid_midpts(self):
        '''
        will return us the Lorentz factor grid midpoints \gamma_{j \pm 1/2}
        it's calculated on gamma_grid_exts o we have \gamma_{-1/2} ... \gamma_{N+1/2}
        with the trick of using the extended grid (\gamma_{-1} and \gamma{N+1}) now
        the grid of delta_gamma  has the same dimension of the (\gamma_{0} and \gamma{N}) grid.
        '''
        return (self.gamma_grid_ext[1:]*self.gamma_grid_ext[:-1])**0.5


    @property
    def delta_gamma(self):
        '''
        will return us the \Delta \gamma as defined after Eq.(5)
        '''
        return self.gamma_grid_midpts[1:] -  self.gamma_grid_midpts[:-1]

    @property
    def nu_synchro_grid(self):
        '''
        return an array of synchrotron frequency corresponing to the gamma_grid
        '''
        return self.gamma_grid**2*(e*self.B)/(2*pi*m_e*c)


    @property
    def N_e_inj(self):
        # particle distribution of injected electrons
        return np.array([self.inj_spectr_norm*gamma**self.inj_spectr_index for gamma in self.gamma_grid])

    @property
    def constant_injection(self):
        '''
        will initialize an array for costant injection
        '''
        return np.array([self.inj_spectr_norm for gamma in self.gamma_grid])

    @property
    def powerlaw_injection(self):
        '''
        will initialize an array for power-law injection
        '''
        return self.N_e_inj
