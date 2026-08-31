import math

from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT


@njit(**njit_dict_no_parallel)
def get_doppler_factor(velocity, mu, enable_full_relativity):
    """
    Calculate the Doppler factor for frame transformation.

    Parameters
    ----------
    velocity : float
        Local packet-frame velocity [cm / s].
    mu : float
        Directional cosine.
    enable_full_relativity : bool
        Whether to use the full relativistic expression.
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = velocity * inv_c
    if not enable_full_relativity:
        return get_doppler_factor_partial_relativity(mu, beta)
    return get_doppler_factor_full_relativity(mu, beta)


@njit(**njit_dict_no_parallel)
def get_doppler_factor_partial_relativity(mu, beta):
    return 1.0 - mu * beta


@njit(**njit_dict_no_parallel)
def get_doppler_factor_full_relativity(mu, beta):
    return (1.0 - mu * beta) / math.sqrt(1 - beta * beta)


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor(velocity, mu, enable_full_relativity):
    """
    Calculate inverse Doppler factor for frame transformation.

    Parameters
    ----------
    velocity : float
        Local packet-frame velocity [cm / s].
    mu : float
        Directional cosine.
    enable_full_relativity : bool
        Whether to use the full relativistic expression.
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = velocity * inv_c
    if not enable_full_relativity:
        return get_inverse_doppler_factor_partial_relativity(mu, beta)
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
    return r_packet.energy * doppler_factor


@njit(**njit_dict_no_parallel)
def angle_aberration_CMF_to_LF(r_packet, time_explosion, mu):
    """
    Converts angle aberration from comoving frame to
    laboratory frame.
    """
    ct = C_SPEED_OF_LIGHT * time_explosion
    beta = r_packet.r / (ct)
    return (mu + beta) / (1.0 + beta * mu)


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
