from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.montecarlo_numba.numba_config import (
    SIGMA_THOMSON,
)

from tardis import constants as const

H = const.h.cgs.value
M_E = const.m_e.cgs.value
K_B = const.k_B.cgs.value
E = const.e.esu.value
C = const.c.cgs.value

FF_OPAC_CONST = (
    (2 * np.pi / (3 * M_E * K_B)) ** 0.5 * 4 * E**6 / (3 * M_E * H * C)
)  # See Eq. 6.1.8 in http://personal.psu.edu/rbc3/A534/lec6.pdf


@njit(**njit_dict_no_parallel)
def chi_electron_calculator(numba_plasma, nu, shell):
    """
    Calculate chi for Thomson scattering

    Parameters
    ----------
    numba_plasma : NumbaPlasma
    nu : float
        Comoving frequency of the r-packet.
    shell : int
        Shell of current r_packet

    Returns
    -------
    numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    """
    return numba_plasma.electron_density[shell] * SIGMA_THOMSON


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
def get_current_bound_free_continua(numba_plasma, nu):
    """
    Determine bound-free continua for which absorption is possible.

    Parameters
    ----------
    numba_plasma : NumbaPlasma
    nu : float
        Comoving frequency of the r-packet.

    Returns
    -------
    numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.
    """
    nu_mins = numba_plasma.photo_ion_nu_threshold_mins
    nu_maxs = numba_plasma.photo_ion_nu_threshold_maxs
    current_continua = np.where(np.logical_and(nu >= nu_mins, nu <= nu_maxs))[0]
    return current_continua


@njit(**njit_dict_no_parallel)
def chi_bf_interpolator(numba_plasma, nu, shell):
    """
    Interpolate the bound-free opacity.

    This function interpolates the tabulated bound-free opacities
    and cross-sections to new frequency values `nu`.

    Parameters
    ----------
    numba_plasma : NumbaPlasma
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

    current_continua = get_current_bound_free_continua(numba_plasma, nu)
    chi_bfs = np.zeros(len(current_continua))
    x_sect_bfs = np.zeros(len(current_continua))
    for i, continuum_id in enumerate(current_continua):
        start = numba_plasma.photo_ion_block_references[continuum_id]
        end = numba_plasma.photo_ion_block_references[continuum_id + 1]
        phot_nus_continuum = numba_plasma.phot_nus[start:end]
        nu_idx = np.searchsorted(phot_nus_continuum, nu)
        interval = phot_nus_continuum[nu_idx] - phot_nus_continuum[nu_idx - 1]
        high_weight = nu - phot_nus_continuum[nu_idx - 1]
        low_weight = phot_nus_continuum[nu_idx] - nu
        chi_bfs_continuum = numba_plasma.chi_bf[start:end, shell]
        chi_bfs[i] = (
            chi_bfs_continuum[nu_idx] * high_weight
            + chi_bfs_continuum[nu_idx - 1] * low_weight
        ) / interval
        x_sect_bfs_continuum = numba_plasma.x_sect[start:end]
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
def chi_ff_calculator(numba_plasma, nu, shell):
    """
    Attributes
    ----------
    numba_plasma : NumbaPlasma
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
        * numba_plasma.ff_opacity_factor[shell]
        / nu**3
        * (1 - np.exp(-H * nu / (K_B * numba_plasma.t_electrons[shell])))
    )
    return chi_ff


@njit(**njit_dict_no_parallel)
def chi_continuum_calculator(numba_plasma, nu, shell):
    """
    Attributes
    ----------
    numba_plasma : NumbaPlasma
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
    ) = chi_bf_interpolator(numba_plasma, nu, shell)
    chi_ff = chi_ff_calculator(numba_plasma, nu, shell)
    return (
        chi_bf_tot,
        chi_bf_contributions,
        current_continua,
        x_sect_bfs,
        chi_ff,
    )
