from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.montecarlo_numba.numba_config import (
    SIGMA_THOMSON,
)


@njit(**njit_dict_no_parallel)
def calculate_tau_electron(electron_density, distance):
    """
    Calculate tau for Thomson scattering

    Parameters
    ----------
    electron_density : float
    distance : float
    """
    return electron_density * SIGMA_THOMSON * distance
