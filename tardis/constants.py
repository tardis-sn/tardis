# astrophysical constants in cgs units (deprecated - suppose to use astropy constants)
import numpy as np

#conversions
erg2ev = 6.24150974e11

#speed of light
c = 2.99792458e10 #cm/s
inverse_c = 1 / c
#masses
mh = 1.67352e-24 #  of hydrogen atom
me = 9.10938188e-28 # electron mass in grams

#luminosities
log_lsun = 33.585009279902458


#charges
e = 4.80320425e-10

#Planck constant
h = 6.6260755e-27

#boltzman constant
kb = 1.38e-16 # erg K^(-1)
kbinev = kb * erg2ev

#stefan-boltzmann constant
sigma_sb = 5.67051e-5

#thomson cross section
sigma_thomson = 6.652486e-25 #cm^(-2)
inverse_sigma_thomson = 1 / sigma_thomson


#sobolev coeff
sobolev_coeff = (np.pi * e ** 2) / (me * c) #(pi*e^2)/(m_e*c)

#
trad_estimator_constant = 0.260944706 * h / kb # (pi**4 / 15)/ (24*zeta(5))

w_estimator_constant = (c ** 2 / (2 * h)) * (15 / np.pi ** 4) * (h / kb) ** 4 / (4 * np.pi)

#atomic mass unit
u = 1.66053886e-24 # atomic mass unit in g
