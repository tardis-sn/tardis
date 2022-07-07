import math

from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

import tardis.montecarlo.montecarlo_numba.numba_config as nc
from tardis.montecarlo.montecarlo_numba.numba_config import (
    C_SPEED_OF_LIGHT,
    MISS_DISTANCE,
    SIGMA_THOMSON,
    CLOSE_LINE_THRESHOLD,
)

from tardis.montecarlo.montecarlo_numba.utils import MonteCarloException
from tardis.montecarlo.montecarlo_numba.r_packet import (
    print_r_packet_properties,
)


def calculate_distance_boundary_2d_cell_spherical(r, mu, r_inner, r_outer):l
    delta_shell = 0
    if mu > 0.0:
        # direction outward
        distance = math.sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (
            r * mu
        )
        delta_shell = 1
    else:
    # going inward
        check = r_inner * r_inner + (r * r * (mu * mu - 1.0))

        if check >= 0.0:
        # hit inner boundary
            distance = -r * mu - math.sqrt(check)
            delta_shell = -1
        else:
            # miss inner boundary
            distance = math.sqrt(
                r_outer * r_outer + ((mu * mu - 1.0) * r * r)
            ) - (r * mu)
            delta_shell = 1
    return distance, delta_shell
