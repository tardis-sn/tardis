import numpy as np
from numba import njit

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)


@njit(**njit_dict_no_parallel)
def piecewise_linear_dvdr(r, shell_id, geometry):
    """
    Velocity at radius r and dv/dr of current shell

    Parameters
    ----------
    r_packet: RPacket
    geometry: Geometry

    Returns
    -------
    v: float, current velocity
    frac: float, dv/dr for current shell
    """
    v_inner = geometry.v_inner[shell_id]
    v_outer = geometry.v_outer[shell_id]
    r_inner = geometry.r_inner[shell_id]
    r_outer = geometry.r_outer[shell_id]
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
    -------
    factor = 1.0 / ((1 - mu * mu) * v / r + mu * mu * dvdr)
    """
    v, dvdr = piecewise_linear_dvdr(r_packet.r, r_packet.current_shell_id, geometry)
    r = r_packet.r
    mu = r_packet.mu
    factor = 1.0 / ((1 - mu * mu) * v / r + mu * mu * dvdr)
    return factor


def depressed_quartic(A, B, C, D, E):

    a = -3.0*B**2.0/(8.0*A**2.0) + C/A
    b = B**3.0/(8.0*A**3.0) - B*C/(2.0*A**2.0) + D/A
    c = -3.0*B**4.0/(256.0*A**4.0) + C*B**2.0/(16.0*A**3.0) - B*D/(4.0*A**2.0) + E/A

    p = -a**2.0/12.0 - c
    q = -a**3.0/108.0 + a*c/3.0 - b**2.0/8.0

    # Handle floating point error to prevent small neg numbers in sqrt
    # Want these terms to exactly cancel if their difference is ~1e-15 as small
    # as their value
    q2p3_thresh = 1e-14 * np.max([np.abs(q**2.0/(4.0)), np.abs(p**3.0/27.0)])
    q2p3_term = q**2.0/(4.0) + p**3.0/27.0
    if np.abs(q2p3_term) <= q2p3_thresh:
        q2p3_term = 0.0

    # Certain solutions will have (physical) complex values of r
    # but the complex component can cancel out in y
    r = -q/2 + np.sqrt(q2p3_term, dtype=complex)
    u = r**(1./3.)

    # Handle case where u=0
    if u == 0.0:
        y = -5.0/6.0 * a - q**(1.0/3.0)
    else:
        y = -5.0/6.0 * a + u - p/(3.0*u)
    y = np.real_if_close(y)
    w = np.sqrt(a + 2.0*y)

    # Handle floating point error to prevent small neg numbers in sqrt
    aybw_term = [-(3.0*a + 2.0*y + 2.0*b/w), -(3.0*a + 2.0*y - 2.0*b/w)]
    aybw_thresh = 1e-13 * np.max([np.abs(3.0*a), np.abs(2.0*y), np.abs(2.0*b/w)])
    for i, term in enumerate(aybw_term):
        if np.abs(term) <= aybw_thresh:
           aybw_term[i] = 0.0 
        elif term < 0.0:
            aybw_term[i] = np.nan

    x1 = -B/(4.0*A) + 0.5*( w + np.sqrt(aybw_term[0]))
    x2 = -B/(4.0*A) + 0.5*( w - np.sqrt(aybw_term[0]))
    x3 = -B/(4.0*A) + 0.5*(-w + np.sqrt(aybw_term[1]))
    x4 = -B/(4.0*A) + 0.5*(-w - np.sqrt(aybw_term[1]))

    roots = (x1, x2, x3, x4)
    return roots
