import numpy as np
import astropy.units as u

from tardis import constants as const
from tardis.energy_input.util import kappa_calculation, SphericalVector

MASS_SI = 28.085 * u.M_p
MASS_FE = 55.845 * u.M_p
    
def compton_opacity_calculation(electron_number_density, energy):
    
    kappa = kappa_calculation(energy)
    
    sigma_T = const.sigma_T.cgs.value

    a = (1. + 2. * kappa)
    
    sigma_KN = 3. / 4. * sigma_T * \
            ((1. + kappa) / kappa ** 3. * ((2. * kappa * (1. + kappa)) / a - np.log(a)) + \
             1. / (2. * kappa) * np.log(a) - (1. + 3 * kappa) / a ** 2.)
    
    return electron_number_density * sigma_KN
    
def photoabsorption_opacity_calculation(energy, ejecta_density, iron_group_fraction):
    
    Si_opacity = 1.16e-24 * (energy / 100.0) ** -3.13 * ejecta_density / MASS_SI.value * \
        (1. - iron_group_fraction)
    
    Fe_opacity = 25.7e-24 * (energy / 100.0) ** -3.0 * ejecta_density / MASS_FE.value * \
        (1. - iron_group_fraction)
    
    return Si_opacity + Fe_opacity

def pair_creation_opacity_calculation(energy, ejecta_density, iron_group_fraction):
    
    Z_Si = 14
    Z_Fe = 26

    Si_proton_ratio = Z_Si ** 2. / MASS_SI.value
    Fe_proton_ratio = Z_Fe ** 2. / MASS_FE.value
    
    multiplier = ejecta_density * (Si_proton_ratio * (1. - iron_group_fraction) + 
        Fe_proton_ratio * iron_group_fraction)
    
    if energy > 1.022e6 and energy < 1.5e6:
        opacity = multiplier * 1.0063 * (energy / 1.0e3 - 1.022) * 1.0e-27
    else:
        opacity = multiplier * (0.0481 + 0.301 * (energy / 1.0e3 - 1.5)) * 1.0e-27
        
    if opacity < 0.:
        opacity = 0.
        
    return opacity