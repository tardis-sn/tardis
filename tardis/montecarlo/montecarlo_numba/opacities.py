import numpy as np
from enum import IntEnum
from numba import int64, float64, boolean
from numba import njit
from numba.experimental import jitclass

import math
from tardis.montecarlo.montecarlo_numba import (
    njit_dict,
    numba_config,
    njit_dict_no_parallel,
)
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.montecarlo.montecarlo_numba.numba_config import (
    CLOSE_LINE_THRESHOLD,
    C_SPEED_OF_LIGHT,
    MISS_DISTANCE,
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
    return electron_density * numba_config.SIGMA_THOMSON * distance