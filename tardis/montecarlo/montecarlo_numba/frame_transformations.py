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
def get_doppler_factor(r, mu, time_explosion):
    inv_c = 1 / C_SPEED_OF_LIGHT
    inv_t = 1 / time_explosion
    beta = r * inv_t * inv_c
    if not numba_config.ENABLE_FULL_RELATIVITY:
        return get_doppler_factor_partial_relativity(mu, beta)
    else:
        return get_doppler_factor_full_relativity(mu, beta)


@njit(**njit_dict_no_parallel)
def get_doppler_factor_partial_relativity(mu, beta):
    return 1.0 - mu * beta


@njit(**njit_dict_no_parallel)
def get_doppler_factor_full_relativity(mu, beta):
    return (1.0 - mu * beta) / math.sqrt(1 - beta * beta)


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor(r, mu, time_explosion):
    """
    Calculate doppler factor for frame transformation

    Parameters
    ----------
    r : float
    mu : float
    time_explosion : float
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    inv_t = 1 / time_explosion
    beta = r * inv_t * inv_c
    if not numba_config.ENABLE_FULL_RELATIVITY:
        return get_inverse_doppler_factor_partial_relativity(mu, beta)
    else:
        return get_inverse_doppler_factor_full_relativity(mu, beta)


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor_partial_relativity(mu, beta):
    return 1.0 / (1.0 - mu * beta)


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor_full_relativity(mu, beta):
    return (1.0 + mu * beta) / math.sqrt(1 - beta * beta)


@njit(**njit_dict_no_parallel)
def calc_packet_energy_full_relativity(r_packet):
    # accurate to 1 / gamma - according to C. Vogl
    return r_packet.energy


@njit(**njit_dict_no_parallel)
def calc_packet_energy(r_packet, distance_trace, time_explosion):
    doppler_factor = 1.0 - (
        (distance_trace + r_packet.mu * r_packet.r)
        / (time_explosion * C_SPEED_OF_LIGHT)
    )
    energy = r_packet.energy * doppler_factor
    return energy


@njit(**njit_dict_no_parallel)
def angle_aberration_CMF_to_LF(r_packet, time_explosion, mu):
    """
    Converts angle aberration from comoving frame to
    laboratory frame.
    """
    ct = C_SPEED_OF_LIGHT * time_explosion
    beta = r_packet.r / (ct)
    return (r_packet.mu + beta) / (1.0 + beta * mu)


@njit(**njit_dict_no_parallel)
def angle_aberration_LF_to_CMF(r_packet, time_explosion, mu):
    """

    c code:
    double beta = rpacket_get_r (packet) * storage->inverse_time_explosion * INVERSE_C;
    return (mu - beta) / (1.0 - beta * mu);
    """
    ct = C_SPEED_OF_LIGHT * time_explosion
    beta = r_packet.r / (ct)
    return (mu - beta) / (1.0 - beta * mu)