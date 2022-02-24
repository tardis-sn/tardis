from numba import njit

from tardis.montecarlo.transport import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.transport.numba_config import (
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
