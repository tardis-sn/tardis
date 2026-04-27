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

@njit(**njit_dict_no_parallel)
def depressed_quartic(A, B, C, D, E):

    a = -3.0*B**2.0/(8.0*A**2.0) + C/A
    b = B**3.0/(8.0*A**3.0) - B*C/(2.0*A**2.0) + D/A
    c = -3.0*B**4.0/(256.0*A**4.0) + C*B**2.0/(16.0*A**3.0) - B*D/(4.0*A**2.0) + E/A

    p = -a**2.0/12.0 - c
    q = -a**3.0/108.0 + a*c/3.0 - b**2.0/8.0

    # Handle floating point error to prevent small neg numbers in sqrt
    # Want these terms to exactly cancel if their difference is ~1e-15 as small
    # as their value
    #q2p3_thresh = 1e-14 * np.max([np.abs(q**2.0/(4.0)), np.abs(p**3.0/27.0)])
    q2p3_term = q**2.0/(4.0) + p**3.0/27.0
    #if np.abs(q2p3_term) <= q2p3_thresh:
    #    q2p3_term = 0.0

    # Certain solutions will have (physical) complex values of r
    # but the complex component can cancel out in y
    q2p3_complex = q2p3_term + 0j
    r = -q/2 + np.sqrt(q2p3_complex)
    u = r**(1./3.)

    # Handle case where u=0
    if u == 0.0:
        y = -5.0/6.0 * a - q**(1.0/3.0)
    else:
        y = -5.0/6.0 * a + u - p/(3.0*u)

    y_real = np.real(y)
    y_imag = np.imag(y)
    # TODO: potentially handle imaginary parts here, but they should be close to zero

    y = y_real
    w = np.sqrt(a + 2.0*y)

    # Handle floating point error to prevent small neg numbers in sqrt
    # connor-mcclellan: note - we lose a lot more precision here than I expect.
    # Edge cases require a threshold as large as 1e-7 to recover the same distance
    # to line as homologous expansion, but I would have expected ~1e-15 would work
    # Somewhat clumsy syntax here is for numba compatibility
    aybw_term, aybw_term_clamped = np.empty(2), np.empty(2)
    aybw_term[0] = -(3.0*a + 2.0*y + 2.0*b/w)
    aybw_term[1] = -(3.0*a + 2.0*y - 2.0*b/w)
    aybw_term_clamped[0] = aybw_term[0]
    aybw_term_clamped[1] = aybw_term[1]
    thresh_factor = 3e-7
    aybw_thresh = thresh_factor * max(np.abs(3.0*a), np.abs(2.0*y), np.abs(2.0*b/w))
    for i in range(2):
        term = aybw_term_clamped[i]
        if np.abs(term) <= aybw_thresh:
            aybw_term_clamped[i] = 0.0
        elif term < 0.0:
            aybw_term_clamped[i] = np.nan

    # Handle no real roots - possibly due to floating point error
    if np.isnan(aybw_term_clamped[0]) and np.isnan(aybw_term_clamped[1]):
        if np.abs(aybw_term[0]) <= np.abs(aybw_term[1]):
            minroot_ind = 0
        else:
            minroot_ind = 1
        thresh_match = np.abs(aybw_term[minroot_ind])/aybw_thresh*thresh_factor
        aybw_term_clamped[minroot_ind] = 0.0
        print(
            "No real roots found in depressed_quartic solver. Smallest term is",
            aybw_term[minroot_ind],
            "(greater than relative threshold for zeroing, which is",
            aybw_thresh,
            "). Bending threshold for this root and forcing to zero.",
            "Threshold is currently", thresh_factor,
            "and would need to be greater than", thresh_match, "to catch this case."
        )

    x1 = -B/(4.0*A) + 0.5*( w + np.sqrt(aybw_term_clamped[0]))
    x2 = -B/(4.0*A) + 0.5*( w - np.sqrt(aybw_term_clamped[0]))
    x3 = -B/(4.0*A) + 0.5*(-w + np.sqrt(aybw_term_clamped[1]))
    x4 = -B/(4.0*A) + 0.5*(-w - np.sqrt(aybw_term_clamped[1]))

    roots = (x1, x2, x3, x4)
    return roots
