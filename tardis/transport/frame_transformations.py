import math

from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT


# TODO: In a future PR collapse to just one function
#       get_dopper_factor(v, mu, enable_full_relativity)
#       and calculate v elsewhere either with v = r / t_explosion
#       or v, _ = piecewise_linear_dvdr(r, geometry)
@njit(**njit_dict_no_parallel)
def get_doppler_factor(r, mu, time_explosion, enable_full_relativity):
    inv_c = 1 / C_SPEED_OF_LIGHT
    inv_t = 1 / time_explosion
    beta = r * inv_t * inv_c
    if not enable_full_relativity:
        return get_doppler_factor_partial_relativity(mu, beta)
    return get_doppler_factor_full_relativity(mu, beta)


@njit(**njit_dict_no_parallel)
def get_doppler_factor_nonhomologous(v, mu, enable_full_relativity):
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = v * inv_c
    if not enable_full_relativity:
        return get_doppler_factor_partial_relativity(mu, beta)
    raise NotImplementedError("Full relativity not implemented for non-homologous mode.")


@njit(**njit_dict_no_parallel)
def get_doppler_factor_from_velocity(v, mu, enable_full_relativity):
    """
    Calculate a Doppler factor from local velocity and propagation angle.

    Parameters
    ----------
    v : float
        Local ejecta velocity.
    mu : float
        Packet propagation angle cosine in the lab frame.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = v * inv_c
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
def get_inverse_doppler_factor(r, mu, time_explosion, enable_full_relativity):
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
    if not enable_full_relativity:
        return get_inverse_doppler_factor_partial_relativity(mu, beta)
    return get_inverse_doppler_factor_full_relativity(mu, beta)


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor_nonhomologous(v, mu, enable_full_relativity):
    """
    Calculate doppler factor for frame transformation

    Parameters
    ----------
    v : float
    mu : float
    enable_full_relativity : bool
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = v * inv_c
    if not enable_full_relativity:
        return get_inverse_doppler_factor_partial_relativity(mu, beta)
    raise NotImplementedError("Full relativity not implemented for non-homologous mode.")


@njit(**njit_dict_no_parallel)
def get_inverse_doppler_factor_from_velocity(v, mu, enable_full_relativity):
    """
    Calculate an inverse Doppler factor from local velocity and angle.

    Parameters
    ----------
    v : float
        Local ejecta velocity.
    mu : float
        Packet propagation angle cosine in the comoving frame.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.
    """
    inv_c = 1 / C_SPEED_OF_LIGHT
    beta = v * inv_c
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


@njit(**njit_dict_no_parallel)
def angle_aberration_CMF_to_LF_from_velocity(v, mu):
    """
    Convert angle aberration from comoving frame to lab frame.

    Parameters
    ----------
    v : float
        Local ejecta velocity.
    mu : float
        Packet propagation angle cosine in the comoving frame.
    """
    beta = v / C_SPEED_OF_LIGHT
    return (mu + beta) / (1.0 + beta * mu)


@njit(**njit_dict_no_parallel)
def angle_aberration_LF_to_CMF_from_velocity(v, mu):
    """
    Convert angle aberration from lab frame to comoving frame.

    Parameters
    ----------
    v : float
        Local ejecta velocity.
    mu : float
        Packet propagation angle cosine in the lab frame.
    """
    beta = v / C_SPEED_OF_LIGHT
    return (mu - beta) / (1.0 - beta * mu)
