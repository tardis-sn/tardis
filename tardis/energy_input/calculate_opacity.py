import numpy as np
import astropy.units as u

from tardis import constants as const
from tardis.energy_input.util import kappa_calculation, SphericalVector

MASS_SI = 28.085 * const.m_p.to(u.g).value
MASS_FE = 55.845 * const.m_p.to(u.g).value
M_P = const.m_p.to(u.g).value
SIGMA_T = const.sigma_T.cgs.value


def compton_opacity_calculation(ejecta_density, energy):

    kappa = kappa_calculation(energy)

    a = 1.0 + 2.0 * kappa

    sigma_KN = (
        3.0
        / 4.0
        * SIGMA_T
        * (
            (1.0 + kappa)
            / kappa ** 3.0
            * ((2.0 * kappa * (1.0 + kappa)) / a - np.log(a))
            + 1.0 / (2.0 * kappa) * np.log(a)
            - (1.0 + 3 * kappa) / a ** 2.0
        )
    )

    return ejecta_density / (M_P * 2) * sigma_KN


def photoabsorption_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):

    Si_opacity = (
        1.16e-24
        * (energy / 100.0) ** -3.13
        * ejecta_density
        / MASS_SI
        * (1.0 - iron_group_fraction)
    )

    Fe_opacity = (
        25.7e-24
        * (energy / 100.0) ** -3.0
        * ejecta_density
        / MASS_FE
        * (1.0 - iron_group_fraction)
    )

    return Si_opacity + Fe_opacity


def pair_creation_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):

    Z_Si = 14
    Z_Fe = 26

    Si_proton_ratio = Z_Si ** 2.0 / MASS_SI
    Fe_proton_ratio = Z_Fe ** 2.0 / MASS_FE

    multiplier = ejecta_density * (
        Si_proton_ratio * (1.0 - iron_group_fraction)
        + Fe_proton_ratio * iron_group_fraction
    )

    if energy > 1.022e3 and energy < 1.5e3:
        opacity = multiplier * 1.0063 * (energy / 1.0e3 - 1.022) * 1.0e-27
    else:
        opacity = (
            multiplier * (0.0481 + 0.301 * (energy / 1.0e3 - 1.5)) * 1.0e-27
        )

    if opacity < 0.0:
        opacity = 0.0

    return opacity
