import math

from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    CLOSE_LINE_THRESHOLD,
    MISS_DISTANCE,
    SIGMA_THOMSON,
)
from tardis.transport.montecarlo.nonhomologous_grid import (
    depressed_quartic,
    piecewise_linear_dvdr,
)
from tardis.transport.montecarlo.utils import MonteCarloException


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
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
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
def calculate_distance_line_nonhomologous(
    rpacket,
    geometry,
    nu_line
):
    """
    Calculate distance until RPacket is in resonance with the next line

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.r_packet.RPacket
    geometry : NumbaRadial1DGeometry
    nu_line : float
        line to check the distance to

    Returns
    -------
    distance (cm)
    """
    r = rpacket.r
    v, _ = piecewise_linear_dvdr(rpacket.r, rpacket.current_shell_id, geometry)
    nu_rest = rpacket.nu
    mu = rpacket.mu

    #TODO: unit check / handling here?
    r_inner = geometry.r_inner[rpacket.current_shell_id]
    r_outer = geometry.r_outer[rpacket.current_shell_id]
    v_inner = geometry.v_inner[rpacket.current_shell_id]
    v_outer = geometry.v_outer[rpacket.current_shell_id]

    # Define useful variables to simplify coefficients
    n = C_SPEED_OF_LIGHT * (1 - nu_line / nu_rest)
    m = (v_outer - v_inner)/(r_outer - r_inner)
    p = 1.0 - mu*mu
    q = v_outer - m*r_outer

    # Characteristic scales for non-dimensionalization
    r0 = r_outer - r_inner
    v0 = v_outer - v_inner

    # Dimensionless quantities to use in the quartic solver - improves floating point accuracy
    rd = r/r0
    nd = n/v0
    md = 1.0
    qd = q/v0

    rd2 = rd*rd
    nd2 = nd*nd
    qd2 = qd*qd

    # Define coefficients of the quartic polynomial
    a = 1.0
    b = -2.0 * nd
    c = nd2 + rd2 * p - qd2
    d = -2.0 * nd * rd2 * p
    e = nd2 * rd2 * p

    # m is the velocity gradient
    # n is the relative line velocity
    # If m and n have the same sign, a doppler shift *may* reach the line in this cell
    # If m and n have opposite signs, the velocity in this cell *cannot* shift the packet towards the line
    beta = v/C_SPEED_OF_LIGHT
    doppler_factor = 1.0 - mu * beta
    comov_nu = nu_rest * doppler_factor

    if (comov_nu - nu_line > 1e-14*nu_line and m > 0.0) or (comov_nu - nu_line < -1e-14*nu_line and m < 0.0):
        # Obtain roots of the quartic polynomial for x (= d_line + r_i \mu_i)
        x = depressed_quartic(a, b, c, d, e)
        # Convert each root x_i to a candidate distance: d = r0*x_i - r*mu
        # Select the smallest positive, finite distance among all four roots.
        distance = MISS_DISTANCE
        for xi in x:
            if math.isnan(xi):
                continue
            d_candidate = r0 * xi - r * mu
            if d_candidate > 0.0 and d_candidate < distance:
                distance = d_candidate
    else:
        distance = MISS_DISTANCE

    if distance <= 0.0 or math.isnan(distance) or distance == MISS_DISTANCE:
        distance = MISS_DISTANCE

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
