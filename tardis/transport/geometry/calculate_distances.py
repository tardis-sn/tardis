import math

from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

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


@njit(**njit_dict_no_parallel)
def calculate_distance_boundary(r, mu, r_inner, r_outer):
    """
    Calculate distance to shell boundary in cm.

    Parameters
    ----------
    r : float
       radial coordinate of the RPacket
    mu : float
       cosine of the direction of movement
    r_inner : float
       inner radius of current shell
    r_outer : float
       outer radius of current shell
    """

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


@njit(**njit_dict_no_parallel)
def calculate_distance_line(
    r_packet,
    comov_nu,
    is_last_line,
    nu_line,
    time_explosion,
    enable_full_relativity,
):
    """
    Calculate distance until RPacket is in resonance with the next line

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    comov_nu : float
        comoving frequency at the CURRENT position of the RPacket
    is_last_line : bool
        return MISS_DISTANCE if at the end of the line list
    nu_line : float
        line to check the distance to
    time_explosion : float
        time since explosion in seconds

    Returns
    -------
    """

    nu = r_packet.nu

    if is_last_line:
        return MISS_DISTANCE

    nu_diff = comov_nu - nu_line

    # for numerical reasons, if line is too close, we set the distance to 0.
    if abs(nu_diff / nu) < CLOSE_LINE_THRESHOLD:
        nu_diff = 0.0

    if nu_diff >= 0:
        distance = (nu_diff / nu) * C_SPEED_OF_LIGHT * time_explosion
    else:
        raise MonteCarloException("nu difference is less than 0.0")

    if enable_full_relativity:
        return calculate_distance_line_full_relativity(
            nu_line, nu, time_explosion, r_packet
        )
    return distance


@njit(**njit_dict_no_parallel)
def calculate_distance_line_full_relativity(
    nu_line, nu, time_explosion, r_packet
):
    # distance = - mu * r + (ct - nu_r * nu_r * sqrt(ct * ct - (1 + r * r * (1 - mu * mu) * (1 + pow(nu_r, -2))))) / (1 + nu_r * nu_r);
    nu_r = nu_line / nu
    ct = C_SPEED_OF_LIGHT * time_explosion
    distance = -r_packet.mu * r_packet.r + (
        ct
        - nu_r
        * nu_r
        * math.sqrt(
            ct * ct
            - (
                1
                + r_packet.r
                * r_packet.r
                * (1 - r_packet.mu * r_packet.mu)
                * (1 + 1.0 / (nu_r * nu_r))
            )
        )
    ) / (1 + nu_r * nu_r)
    return distance


@njit(**njit_dict_no_parallel)
def calculate_distance_electron(electron_density, tau_event):
    """
    Calculate distance to Thomson Scattering

    Parameters
    ----------
    electron_density : float
    tau_event : float
    """
    # add full_relativity here
    return tau_event / (electron_density * SIGMA_THOMSON)
