import numpy as np
import radioactivedecay as rd
from numba import njit

from tardis import constants as const
from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)
from tardis.transport.montecarlo.numba_config import (
    SIGMA_THOMSON,
)

H = const.h.cgs.value
M_E = const.m_e.cgs.value
K_B = const.k_B.cgs.value
E = const.e.esu.value
C = const.c.cgs.value
M_P = const.m_p.cgs.value
MASS_SI = rd.Nuclide("Si-28").atomic_mass * M_P
MASS_FE = rd.Nuclide("Fe-56").atomic_mass * M_P
SIGMA_T = const.sigma_T.cgs.value
FINE_STRUCTURE = const.alpha.value
ELECTRON_MASS_ENERGY_KEV = (const.m_e * const.c**2.0).to("keV").value

FF_OPAC_CONST = (
    (2 * np.pi / (3 * M_E * K_B)) ** 0.5 * 4 * E**6 / (3 * M_E * H * C)
)  # See Eq. 6.1.8 in http://personal.psu.edu/rbc3/A534/lec6.pdf


@njit(**njit_dict_no_parallel)
def kappa_calculation(energy):
    """
    Calculates kappa for various other calculations
    i.e. energy normalized to electron rest energy
    511.0 KeV

    Parameters
    ----------
    energy : float

    Returns
    -------
    kappa : float

    """
    return energy / ELECTRON_MASS_ENERGY_KEV


@njit(**njit_dict_no_parallel)
def chi_electron_calculator(opacity_state, nu, shell):
    """
    Calculate chi for Thomson scattering

    Parameters
    ----------
    opacity_state : OpacityState
    nu : float
        Comoving frequency of the r-packet.
    shell : int
        Shell of current r_packet

    Returns
    -------
    numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    """
    return opacity_state.electron_density[shell] * SIGMA_THOMSON


@njit(**njit_dict_no_parallel)
def calculate_tau_electron(electron_density, distance):
    """
    Calculate tau for Thomson scattering

    Parameters
    ----------
    electron_density : float
    distance : float

    Returns
    -------
    tau_electron : float
        tau for thomson scattering
    """
    return electron_density * SIGMA_THOMSON * distance


@njit(**njit_dict_no_parallel)
def get_current_bound_free_continua(opacity_state, nu):
    """
    Determine bound-free continua for which absorption is possible.

    Parameters
    ----------
    opacity_state : OpacityState
    nu : float
        Comoving frequency of the r-packet.

    Returns
    -------
    numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    """
    nu_mins = opacity_state.photo_ion_nu_threshold_mins
    nu_maxs = opacity_state.photo_ion_nu_threshold_maxs
    current_continua = np.where(np.logical_and(nu >= nu_mins, nu <= nu_maxs))[0]
    return current_continua


@njit(**njit_dict_no_parallel)
def chi_bf_interpolator(opacity_state, nu, shell):
    """
    Interpolate the bound-free opacity.

    This function interpolates the tabulated bound-free opacities
    and cross-sections to new frequency values `nu`.

    Parameters
    ----------
    opacity_state : OpacityState
    nu : float, dtype float
        Comoving frequency of the r-packet.
    shell : int, dtype float
        Current computational shell.

    Returns
    -------
    chi_bf_tot : float
        Total bound-free opacity at frequency `nu`.
    chi_bf_contributions : numpy.ndarray, dtype float
        Cumulative distribution function of the contributions of the
        individual bound free continua to the total bound-free opacity.
    current_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    x_sect_bfs : numpy.ndarray, dtype float
        Photoionization cross-sections of all bound-free continua for
        which absorption is possible for frequency `nu`.
    """
    current_continua = get_current_bound_free_continua(opacity_state, nu)
    chi_bfs = np.zeros(len(current_continua))
    x_sect_bfs = np.zeros(len(current_continua))
    for i, continuum_id in enumerate(current_continua):
        start = opacity_state.photo_ion_block_references[continuum_id]
        end = opacity_state.photo_ion_block_references[continuum_id + 1]
        phot_nus_continuum = opacity_state.phot_nus[start:end]
        nu_idx = np.searchsorted(phot_nus_continuum, nu)
        interval = phot_nus_continuum[nu_idx] - phot_nus_continuum[nu_idx - 1]
        high_weight = nu - phot_nus_continuum[nu_idx - 1]
        low_weight = phot_nus_continuum[nu_idx] - nu
        chi_bfs_continuum = opacity_state.chi_bf[start:end, shell]
        chi_bfs[i] = (
            chi_bfs_continuum[nu_idx] * high_weight
            + chi_bfs_continuum[nu_idx - 1] * low_weight
        ) / interval
        x_sect_bfs_continuum = opacity_state.x_sect[start:end]
        x_sect_bfs[i] = (
            x_sect_bfs_continuum[nu_idx] * high_weight
            + x_sect_bfs_continuum[nu_idx - 1] * low_weight
        ) / interval

    chi_bf_contributions = chi_bfs.cumsum()

    # If we are outside the range of frequencies
    # for which we have photo-ionization cross sections
    # we will have no local continuua and therefore
    # no bound-free interactions can occur
    # so we set the bound free opacity to zero
    if len(current_continua) == 0:
        chi_bf_tot = 0.0
    else:
        chi_bf_tot = chi_bf_contributions[-1]
        chi_bf_contributions /= chi_bf_tot

    return (
        chi_bf_tot,
        chi_bf_contributions,
        current_continua,
        x_sect_bfs,
    )


@njit(**njit_dict_no_parallel)
def chi_ff_calculator(opacity_state, nu, shell):
    """
    Attributes
    ----------
    opacity_state : OpacityState
    nu : float64
        Comoving frequency of the r_packet
    shell : int64
        Current shell id of the r_packet

    Returns
    -------
        chi_ff : float64
            Free Free opacity
    """
    chi_ff = (
        FF_OPAC_CONST
        * opacity_state.ff_opacity_factor[shell]
        / nu**3
        * (1 - np.exp(-H * nu / (K_B * opacity_state.t_electrons[shell])))
    )
    return chi_ff


@njit(**njit_dict_no_parallel)
def chi_continuum_calculator(opacity_state, nu, shell):
    """
    Attributes
    ----------
    opacity_state : OpacityState
    nu : float64
        Comoving frequency of the r_packet
    shell : int64
        Current shell id of the r_packet

    Returns
    -------
    chi_bf_tot : float
        Total bound-free opacity at frequency `nu`.
    chi_bf_contributions : numpy.ndarray, dtype float
        Cumulative distribution function of the contributions of the
        individual bound free continua to the total bound-free opacity.
    current_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    x_sect_bfs : numpy.ndarray, dtype float
        Photoionization cross-sections of all bound-free continua for
        which absorption is possible for frequency `nu`.
    chi_ff : float64
            Free Free opacity at frequency `nu`
    """
    (
        chi_bf_tot,
        chi_bf_contributions,
        current_continua,
        x_sect_bfs,
    ) = chi_bf_interpolator(opacity_state, nu, shell)
    chi_ff = chi_ff_calculator(opacity_state, nu, shell)
    return (
        chi_bf_tot,
        chi_bf_contributions,
        current_continua,
        x_sect_bfs,
        chi_ff,
    )


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
    \\rho / 2 p_m \\times 3/4 \\sigma_T ((1 + \\kappa) / \\kappa^3
    ((2\\kappa(1 + \\kappa)) / (1 + 2\\kappa) - \\ln(1 + 2\\kappa) + 1/(2\\kappa) \\ln(1 + 2\\kappa)
    - (1 + 3\\kappa) / (1 + 2\\kappa)^2)
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
