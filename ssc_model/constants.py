# variables with class-wise scope
# remember to use
# from constants import *
# and don't use any of the following symbles, this is to simplify formulas

from astropy import constants as const

sigma_T =  const.sigma_T.to('cm2').value
c = const.c.to('cm/s').value
m_e = const.m_e.to('g').value
E_rest = (const.m_e * const.c**2).to('eV').value
e = const.e.esu.value
h = const.h.to('erg*s').value
# critical magnetic field for the analytical solution
B_c = 4.414*1e13 # in Gauss
