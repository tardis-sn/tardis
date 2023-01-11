import numpy as np
import astropy.units as u
import radioactivedecay as rd
from numba import njit

from tardis import constants as const
from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.energy_input.util import kappa_calculation

# Force cgs
M_P = const.m_p.to(u.g).value
MASS_SI = rd.Nuclide("Si-28").atomic_mass * M_P
MASS_FE = rd.Nuclide("Fe-56").atomic_mass * M_P
SIGMA_T = const.sigma_T.cgs.value
FINE_STRUCTURE = const.alpha.value


@njit(**njit_dict_no_parallel)
def compton_opacity_partial(energy, fraction):
    """Partial Compton scattering opacity, from artis file
    gamma.cc

    Parameters
    ----------
    energy : float
        packet energy in terms of electron rest energy
    fraction : float
        1 + 2 * packet energy

    Returns
    -------
    np.float64
        Compton scattering opacity
    """
    term_1 = (
        ((energy**2.0) - (2.0 * energy) - 2.0)
        * np.log(fraction)
        / energy
        / energy
    )
    term_2 = (((fraction**2.0) - 1.0) / (fraction**2.0)) / 2.0
    term_3 = ((fraction - 1.0) / energy) * (
        (1 / energy) + (2.0 / fraction) + (1.0 / (energy * fraction))
    )
    return 3.0 * SIGMA_T * (term_1 + term_2 + term_3) / (8 * energy)


@njit(**njit_dict_no_parallel)
def compton_opacity_calculation(energy, electron_density):
    """Calculate the Compton scattering opacity for a given energy
    (Rybicki & Lightman, 1979)

    $
    \\rho / 2 p_m \\times 3/4 \\sigma_T ((1 + \kappa) / \kappa^3
    ((2\kappa(1 + \kappa)) / (1 + 2\kappa) - \ln(1 + 2\kappa) + 1/(2\kappa) \ln(1 + 2\kappa)
    - (1 + 3\kappa) / (1 + 2\kappa)^2)
    $

    Parameters
    ----------
    energy : float
        The energy of the photon
    ejecta_density : float
        The density of the ejecta

    Returns
    -------
    float
        The Compton scattering opacity
    """
    kappa = kappa_calculation(energy)

    a = 1.0 + 2.0 * kappa

    sigma_KN = (
        3.0
        / 4.0
        * SIGMA_T
        * (
            (1.0 + kappa)
            / kappa**3.0
            * ((2.0 * kappa * (1.0 + kappa)) / a - np.log(a))
            + 1.0 / (2.0 * kappa) * np.log(a)
            - (1.0 + 3 * kappa) / a**2.0
        )
    )

    return electron_density * sigma_KN


@njit(**njit_dict_no_parallel)
def photoabsorption_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):
    """Calculates photoabsorption opacity for a given energy
    Approximate treatment from Ambwani & Sutherland (1988)
    Magic numbers are from the approximate treatment

    Parameters
    ----------
    energy : float
        Photon energy
    ejecta_density : float
        The density of the ejecta
    iron_group_fraction : float
        Fraction of iron group elements in the shell

    Returns
    -------
    float
        Photoabsorption opacity
    """
    si_opacity = (
        1.16e-24
        * (energy / 100.0) ** -3.13
        * ejecta_density
        / MASS_SI
        * (1.0 - iron_group_fraction)
    )

    fe_opacity = (
        25.7e-24
        * (energy / 100.0) ** -3.0
        * ejecta_density
        / MASS_FE
        * iron_group_fraction
    )

    return si_opacity + fe_opacity


@njit(**njit_dict_no_parallel)
def photoabsorption_opacity_calculation_kasen(
    energy, number_density, proton_count
):
    """Calculates photoabsorption opacity for a given energy
    Approximate treatment from Kasen et al. (2006)

    Parameters
    ----------
    energy : float
        Photon energy
    number_density : float
        The number density of the ejecta for each atom
    proton_count : float
        Number of protons for each atom in the ejecta

    Returns
    -------
    float
        Photoabsorption opacity
    """
    kappa = kappa_calculation(energy)

    opacity = (FINE_STRUCTURE**4.0) * 8.0 * np.sqrt(2) * (kappa**-3.5)
    # Note- this should actually be atom_number_density * (atom_proton_number ** 5)
    return (
        SIGMA_T
        * opacity
        * np.sum((number_density / proton_count) * proton_count**5)
    )


@njit(**njit_dict_no_parallel)
def pair_creation_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):
    """Calculates pair creation opacity for a given energy
    Approximate treatment from Ambwani & Sutherland (1988)

    Parameters
    ----------
    energy : float
        Photon energy
    ejecta_density : float
        The density of the ejecta
    iron_group_fraction : float
        Fraction of iron group elements in the shell

    Returns
    -------
    float
        Pair creation opacity
    """
    z_si = 14
    z_fe = 26

    si_proton_ratio = z_si**2.0 / MASS_SI
    fe_proton_ratio = z_fe**2.0 / MASS_FE

    multiplier = ejecta_density * (
        si_proton_ratio * (1.0 - iron_group_fraction)
        + fe_proton_ratio * iron_group_fraction
    )

    # Conditions prevent divide by zero
    # Ambwani & Sutherland (1988)
    if energy > 1022 and energy < 1500:
        opacity = multiplier * 1.0063 * (energy / 1000 - 1.022) * 1.0e-27
    elif energy >= 1500:
        opacity = (
            multiplier * (0.0481 + 0.301 * (energy / 1000 - 1.5)) * 1.0e-27
        )
    else:
        opacity = 0

    return opacity


@njit(**njit_dict_no_parallel)
def pair_creation_opacity_artis(energy, ejecta_density, iron_group_fraction):
    """Calculates pair creation opacity for a given energy
    Approximate treatment from Ambwani & Sutherland (1988)
    as implemented in ARTIS

    Parameters
    ----------
    energy : float
        Photon energy
    ejecta_density : float
        The density of the ejecta
    iron_group_fraction : float
        Fraction of iron group elements in the shell

    Returns
    -------
    float
        Pair creation opacity
    """
    # Conditions prevent divide by zero
    # Ambwani & Sutherland (1988)
    if energy > 1022:
        if energy > 1500:
            opacity_si = (0.0481 + (0.301 * (energy - 1500))) * 196.0e-27
            opacity_fe = (0.0481 + (0.301 * (energy - 1500))) * 784.0e-27
        else:
            opacity_si = 1.0063 * (energy - 1022) * 196.0e-27
            opacity_fe = 1.0063 * (energy - 1022) * 784.0e-27

        opacity_si *= ejecta_density / M_P / 28
        opacity_fe *= ejecta_density / M_P / 56

        opacity = (opacity_fe * iron_group_fraction) + (
            opacity_si * (1.0 - iron_group_fraction)
        )
    else:
        opacity = 0

    return opacity
