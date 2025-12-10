from collections import namedtuple

import astropy.units as units
import numpy as np

import tardis.constants as const

Constants = namedtuple("Constants", ["c_einstein", "C0_ff", "c0_reg", "I_H"])

c_einstein = (
    4.0 * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)
).value
C0_ff = 1.426e-27  # in cgs units (see Osterbrock 1974)
c_0_regemorter = 5.465e-11
I_H = 13.598433770784 * units.eV.to(units.erg)

continuum_constants = Constants(
    c_einstein=c_einstein, c0_reg=c_0_regemorter, I_H=I_H, C0_ff=C0_ff
)
