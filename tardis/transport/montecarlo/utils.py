import numpy as np
from numba import njit
from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)


class MonteCarloException(ValueError):
    pass


@njit(**njit_dict_no_parallel)
def get_random_mu():
    return 2.0 * np.random.random() - 1.0
