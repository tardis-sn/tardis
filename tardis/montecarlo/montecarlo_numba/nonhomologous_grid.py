import numpy as np
from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)


@njit(**njit_dict_no_parallel)
def velocity(r_packet, numba_model):
    """
    Velocity at radius r

    Parameters
    ----------
    r_packet: RPacket
    numba_model: NumbaModel

    Returns
    -----------
    v: float, current velocity
    """
    shell_id = r_packet.current_shell_id
    v_inner = numba_model.v_inner[shell_id]
    v_outer = numba_model.v_outer[shell_id]
    r_inner = numba_model.r_inner[shell_id]
    r_outer = numba_model.r_outer[shell_id]
    r = r_packet.r
    frac = (v_outer - v_inner) / (r_outer - r_inner)
    return v_inner + frac * (r - r_inner)


@njit(**njit_dict_no_parallel)
def dv_dr(r_packet, numba_model):
    """
    dv/dr for the current shell

    Parameters
    ----------
    r_packet: RPacket
    numba_model: NumbaModel

    Returns
    -----------
    dvdr: float, dv/dr of the current shell
    """
    shell_id = r_packet.current_shell_id
    v_inner = numba_model.v_inner[shell_id]
    v_outer = numba_model.v_outer[shell_id]
    r_inner = numba_model.r_inner[shell_id]
    r_outer = numba_model.r_outer[shell_id]
    return (v_outer - v_inner) / (r_outer - r_inner)


@njit(**njit_dict_no_parallel)
def tau_sobolev_factor(r_packet, numba_model, numba_plasma):
    """
    The angle and velocity dependent Tau Sobolev factor. Is called when ENABLE_NONHOMOLOGOUS_EXPANSION is set to True.

    Parameters
    ----------
    r_packet: RPacket
    numba_model: NumbaModel
    numba_plasma: NumbaPlasma

    Returns
    -----------
    factor = //put the equation here
    """
    shell_id = r_packet.current_shell_id

    v = velocity(r_packet, numba_model)
    r = r_packet.r
    dvdr = dv_dr(r_packet, numba_model)
    mu = r_packet.mu
    factor = 1.0 / ((1 - mu * mu) * v / r + mu * mu * dvdr)
    return factor


@njit(**njit_dict_no_parallel)
def roots(a, b, c, d, e, threshold):
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
    delta = (
        256 * a**3 * e**3
        - 192 * a**2 * b * d * e**2
        - 128 * a**2 * c**2 * e**2
        + 144 * a**2 * c * d**2 * e
        - 27 * a**2 * d**4
        + 144 * a * b**2 * c * e**2
        - 6 * a * b**2 * d**2 * e
        - 80 * a * b * c**2 * d * e
        + 18 * a * b * c * d**3
        + 16 * a * c**4 * e
        - 4 * a * c**3 * d**2
        - 27 * b**4 * e**2
        + 18 * b**3 * c * d * e
        - 4 * b**3 * d**3
        - 4 * b**2 * c**3 * e
        + b**2 * c**2 * d**2
    )
    P = 8 * a * c - 3 * b**2
    R = b**3 + 8 * d * a**2 - 4 * a * b * c
    delta_0 = c**2 - 3 * b * d + 12 * a * e
    delta_1 = (
        2 * c**3
        - 9 * b * c * d
        + 27 * a * d**2
        + 27 * b**2 * e
        - 72 * a * c * e
    )
    D = (
        64 * a**3 * e
        - 16 * a**2 * c**2
        + 16 * a * b**2 * c
        - 16 * a**2 * b * d
        - 3 * b**4
    )
    p = (8 * a * c - 3 * b**2) / (8 * a**2)
    q = (b**3 - 4 * a * b * c + 8 * a**2 * d) / (8 * a**3)
    x_1, x_2, x_3, x_4 = None, None, None, None
    if delta < 0:

        Q = ((delta_1 + np.sqrt(delta_1**2 - 4 * delta_0**3)) / 2) ** (
            1 / 3
        )
        S = np.sqrt(-2 / 3 * p + 1 / (3 * a) * (Q + delta_0 / Q)) / 2

        if -4 * S**2 - 2 * p - q / S >= 0:
            x_1 = (
                -b / (4 * a) + S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_2 = (
                -b / (4 * a) + S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_3 = None
            x_4 = None
        if -4 * S**2 - 2 * p + q / S >= 0:
            x_3 = (
                -b / (4 * a) - S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
            x_4 = (
                -b / (4 * a) - S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
            x_1 = None
            x_2 = None
    elif delta > 0:
        if P < 0 and D < 0:
            phi = np.arccos(delta_1 / (2 * np.sqrt(delta_0**3)))
            S = (
                np.sqrt(
                    -2 / 3 * p
                    + 2 / (3 * a) * np.sqrt(delta_0) * np.cos(phi / 3)
                )
                / 2
            )
            x_1 = (
                -b / (4 * a) + S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_2 = (
                -b / (4 * a) + S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_3 = (
                -b / (4 * a) - S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
            x_4 = (
                -b / (4 * a) - S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
        elif P > 0 and D > 0:
            pass
    else:
        if P < 0 and D < 0 and delta_0 != 0:
            phi = np.arccos(delta_1 / (2 * np.sqrt(delta_0**3)))
            S = (
                np.sqrt(
                    -2 / 3 * p
                    + 2 / (3 * a) * np.sqrt(delta_0) * np.cos(phi / 3)
                )
                / 2
            )
            x_1 = (
                -b / (4 * a) + S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_2 = (
                -b / (4 * a) + S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_3 = (
                -b / (4 * a) - S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
            x_4 = (
                -b / (4 * a) - S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
        elif (D > 0) or (P > 0 and (D != 0 or R != 0)):
            phi = np.arccos(delta_1 / (2 * np.sqrt(delta_0**3)))
            S = (
                np.sqrt(
                    -2 / 3 * p
                    + 2 / (3 * a) * np.sqrt(delta_0) * np.cos(phi / 3)
                )
                / 2
            )
            if -4 * S**2 - 2 * p - q / S >= 0:
                x_1 = (
                    -b / (4 * a)
                    + S
                    + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
                )
                x_2 = (
                    -b / (4 * a)
                    + S
                    - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
                )
            if -4 * S**2 - 2 * p + q / S >= 0:
                x_3 = (
                    -b / (4 * a)
                    - S
                    + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
                )
                x_4 = (
                    -b / (4 * a)
                    - S
                    - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
                )
        elif delta_0 == 0 and D != 0:
            phi = np.arccos(delta_1 / (2 * np.sqrt(delta_0**3)))
            S = (
                np.sqrt(
                    -2 / 3 * p
                    + 2 / (3 * a) * np.sqrt(delta_0) * np.cos(phi / 3)
                )
                / 2
            )
            x_1 = (
                -b / (4 * a) + S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_2 = (
                -b / (4 * a) + S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
            )
            x_3 = (
                -b / (4 * a) - S + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
            x_4 = (
                -b / (4 * a) - S - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
            )
        elif D == 0:
            if P < 0:
                phi = np.arccos(delta_1 / (2 * np.sqrt(delta_0**3)))
                S = (
                    np.sqrt(
                        -2 / 3 * p
                        + 2 / (3 * a) * np.sqrt(delta_0) * np.cos(phi / 3)
                    )
                    / 2
                )
                x_1 = (
                    -b / (4 * a)
                    + S
                    + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
                )
                x_2 = (
                    -b / (4 * a)
                    + S
                    - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p - q / S)
                )
                x_3 = (
                    -b / (4 * a)
                    - S
                    + 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
                )
                x_4 = (
                    -b / (4 * a)
                    - S
                    - 1 / 2 * np.sqrt(-4 * S**2 - 2 * p + q / S)
                )
            elif P > 0 and R == 0:
                pass
            elif delta_0 == 0:
                x_1 = -b / (4 * a)
    roots = []
    if x_1 != None:
        x_1 -= threshold
        if x_1 > 0:
            roots.append(x_1)
    if x_2 != None:
        x_2 -= threshold
        if x_2 > 0:
            roots.append(x_2)
    if x_3 != None:
        x_3 -= threshold
        if x_3 > 0:
            roots.append(x_3)
    if x_4 != None:
        x_4 -= threshold
        if x_4 > 0:
            roots.append(x_4)

    return roots
