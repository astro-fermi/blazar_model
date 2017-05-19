from __future__ import division
from math import pi, sqrt, exp, log
import numpy as np
from scipy import special, integrate
from constants import *
import astropy.units as u
import astropy.constants as const
import naima

class numerics:
    '''
    Numeric procedures needed to evaluate SSC emission.
    class to be intialized with a model object (as defined in module.py)
    Reference is Chiaberge & Ghisellini (1998) MNRAS, 306,551 (1999).
    '''


    def __init__(self, model_object):
        # initialize the model object we will use from now on
        self.model = model_object


    def ChaCoo_tridiag_matrix(self):
        '''
        Implementing tridiagonal matrix of Eq.(9) of Reference
        '''
        # the dimension of our problem is
        N = len(self.model.delta_gamma)
        ChaCoo_matrix = np.zeros((N,N), float)

        # function for the cooling rate
        def cool_rate(U_B, U_rad, gamma):
            '''
            cooling rate from Eq.(2)
            '''
            return 4/3 * sigma_T/(m_e*c) * (U_B + U_rad) * gamma**2

        # loop on the energy to fill the matrix
        for i in range(N):
            gamma_minus_half = self.model.gamma_grid_midpts[i]                              # \gamma_{j-1/2} of Reference
            gamma_plus_half = self.model.gamma_grid_midpts[i] + self.model.delta_gamma[i]   # \gamma_{j+1/2} of Reference
            delta_gamma = self.model.delta_gamma[i]
            delta_t = self.model.delta_t
            t_esc = self.model.t_esc
            U_B = self.model.U_B
            # when we'll introduce synchrotron radiation U_rad will be an array
            # with a value for each point of the gamma_grid (the extremes of integration depend on gamma)
            U_rad = self.model.U_rad[i]
            # Eq.s (10) of Reference
            V2 = 1 + delta_t/t_esc + \
                     (delta_t * cool_rate(U_B, U_rad, gamma_minus_half))/delta_gamma
            V3 = - (delta_t * cool_rate(U_B, U_rad, gamma_plus_half))/delta_gamma
            # let's loop on another dimension to fill the diagonal of the matrix
            for j in range(N):
                if j == i:
                    ChaCoo_matrix[i,j] = V2
                if j == i+1:
                    ChaCoo_matrix[i,j] = V3

        # Chang Cooper boundaries condition, Eq.(18) Park Petrosian, 1996
        # http://adsabs.harvard.edu/full/1996ApJS..103..255P
        ChaCoo_matrix[N-2,N-1] = 0.
        return ChaCoo_matrix


    def synchrotron(self, N_e):
        # we plug now the obtained electron distribution in naima TableModel
        # remultiply by Volume and maximum time of evolution to get N_e differential in Energy
        N_e_differential = N_e * 4/3*pi*self.model.R**3 * self.model.time_max
        electron_density = naima.models.TableModel( self.model.energy_grid * u.eV,
                                                    N_e_differential * u.Unit('1/erg'),
                                                    amplitude = 1)

        SYN = naima.models.Synchrotron(electron_density, B = self.model.B * u.G)

        return SYN


    def inverse_compton(self, N_e):

        N_e_differential = N_e * 4/3*pi*self.model.R**3 * self.model.time_max
        electron_density = naima.models.TableModel( self.model.energy_grid * u.eV,
                                                    N_e_differential * u.Unit('1/erg'),
                                                    amplitude = 1)

        # Define energy array for synchrotron seed photon field and compute
        # Synchroton luminosity by setting distance to 0.
        energy = np.logspace(-7, 15, 100) * u.eV
        Lsy = self.synchrotron(N_e).flux(energy, distance=0*u.cm)

        # Define source radius and compute photon density
        R =  self.model.R * u.cm
        phn_sy = Lsy / (4 * np.pi * R**2 * const.c) * 2.26

        # Create IC instance with CMB and synchrotron seed photon fields:
        IC = naima.models.InverseCompton(electron_density, seed_photon_fields=[['SSC', energy, phn_sy]])

        return IC


    def evolve(self):
        '''
        Evolving injected spectrum solving iteratively Eq.(9).
        We will calculate the synchrotron emissivity with the romberg integration
        and will update the model.U_rad parameter.
        Options only_synchrotron_cooling is for test
        '''
        # injected spectrum
        N_e = self.model.N_e_inj
        # injecton term, to be added each delta_t up to the maximum injection time
        # specified by model.inj_time
        Q_e = self.model.powerlaw_injection
        delta_t = self.model.delta_t

        # time grid loop
        time_past_injection = 0.

        for time in self.model.time_grid:
            # here N^{i+1} of Reference --> N_e_tmp
            # here N^{i} of Reference --> N_e

            # solve the system with Eq.(11):
            if time_past_injection <= self.model.inj_time:
                N_e_tmp = np.linalg.solve(self.ChaCoo_tridiag_matrix(), N_e + Q_e*delta_t)
            # no more injection after the established t_inj
            else:
                N_e_tmp = np.linalg.solve(self.ChaCoo_tridiag_matrix(), N_e)

            # swap!, now we are moving to the following istant of time
            # N^{i+1} --> N^{i} and restart
            N_e = N_e_tmp
            # update the time past injection
            print 't after injection: ', time_past_injection/self.model.crossing_time, ' crossing time'
            time_past_injection += delta_t

        return N_e
