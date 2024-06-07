import numpy as np
from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)


@njit(**njit_dict_no_parallel)
def velocity_dvdr(r_packet, geometry):
    """
    Velocity at radius r and dv/dr of current shell

    Parameters
    ----------
    r_packet: RPacket
    geometry: Geometry

    Returns
    -----------
    v: float, current velocity
    frac: float, dv/dr for current shell
    """
    shell_id = r_packet.current_shell_id
    v_inner = geometry.v_inner[shell_id]
    v_outer = geometry.v_outer[shell_id]
    r_inner = geometry.r_inner[shell_id]
    r_outer = geometry.r_outer[shell_id]
    r = r_packet.r
    frac = (v_outer - v_inner) / (r_outer - r_inner)
    return v_inner + frac * (r - r_inner), frac


@njit(**njit_dict_no_parallel)
def tau_sobolev_factor(r_packet, geometry):
    """
    The angle and velocity dependent Tau Sobolev factor component. Is called when ENABLE_NONHOMOLOGOUS_EXPANSION is set to True.

    Note: to get Tau Sobolev, this needs to be multiplied by tau_sobolevs found from plasma
    Parameters
    ----------
    r_packet: RPacket
    geometry: Geometry

    Returns
    -----------
    factor = 1.0 / ((1 - mu * mu) * v / r + mu * mu * dvdr)
    """

    v, dvdr = velocity_dvdr(r_packet, geometry)
    r = r_packet.r
    mu = r_packet.mu
    factor = 1.0 / ((1 - mu * mu) * v / r + mu * mu * dvdr)
    return factor


# @njit(**njit_dict_no_parallel)
def quartic_roots(a, b, c, d, e, threshold):
    """
    Solves ax^4 + bx^3 + cx^2 + dx + e = 0, for the real roots greater than the threshold returns (x - threshold).
    Uses: https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots

    Parameters
    -----------
    a, b, c, d, e: coefficients of the equations ax^4 + bx^3 + cx^2 + dx + e = 0, float
    threshold: lower needed limit on roots, float
    Returns
    -----------
    roots: real positive roots of ax^4 + bx^3 + cx^2 + dx + e = 0

    """
    roots = np.roots((a, b, c, d, e))
    roots = [root for root in roots if isinstance(root, float)]
    roots = [root for root in roots if root > threshold]

    return roots
