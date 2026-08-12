import math

import numpy as np
import numpy.typing as npt
from numba import njit

from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
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
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket
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
        return 0.0  # Double check that this is valid. Catches full relativity case.

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
    rpacket: RPacket,
    geometry: NumbaNonhomologousRadial1DGeometry,
    nu_line: float,
    minimum_distance: float = 0.0,
    maximum_distance: float = MISS_DISTANCE,
) -> float:
    """
    Calculate distance until RPacket is in resonance with the next line

    Candidate roots must lie after ``minimum_distance`` and no farther than
    ``maximum_distance`` along the current packet trajectory.

    Returns
    -------
    distance (cm)
    """
    #TODO: unit check / handling here?
    r_inner = geometry.r_inner[rpacket.current_shell_id]
    r_outer = geometry.r_outer[rpacket.current_shell_id]
    v_inner = geometry.v_inner[rpacket.current_shell_id]
    v_outer = geometry.v_outer[rpacket.current_shell_id]

    r = rpacket.r
    dvdr = geometry.velocity_gradient[rpacket.current_shell_id]
    nu_rest = rpacket.nu
    mu = rpacket.mu

    # Define useful variables to simplify coefficients
    n = C_SPEED_OF_LIGHT * (1 - nu_line / nu_rest)
    m = dvdr
    p = 1.0 - mu*mu
    q = v_outer - m*r_outer

    # Characteristic scales for non-dimensionalization
    r0 = r_outer - r_inner
    v0 = v_outer - v_inner
    distance_tolerance = CLOSE_LINE_THRESHOLD * r0
    if v0 == 0.0:
        impact_parameter_squared = r * r * p
        if q != 0.0 and n * n < q * q:
            x_squared = (
                n * n * impact_parameter_squared / (q * q - n * n)
            )
            x_root = math.copysign(math.sqrt(x_squared), n / q) / r0
            x = (x_root, math.nan, math.nan, math.nan)
        else:
            x = (math.nan, math.nan, math.nan, math.nan)
    else:
        # Dimensionless quantities to use in the quartic solver - improves floating point accuracy
        rd = r / r0
        nd = n / v0
        md = 1.0  # m/(v0/r0) # dimensionless m will always be 1
        qd = q / v0

        md2 = md * md
        rd2 = rd * rd
        nd2 = nd * nd
        qd2 = qd * qd

        # Define coefficients of the quartic polynomial
        a = md2
        b = -2.0 * nd * md
        c = nd2 + md * rd2 * p - qd2
        d = -2.0 * nd * md * rd2 * p
        e = nd2 * rd2 * p

        # Obtain roots of the quartic polynomial for x (= d_line + r_i \mu_i)
        x = depressed_quartic(a, b, c, d, e)

    # Convert each root x_i to a candidate distance: d = r0*x_i - r*mu
    # Select the nearest root that satisfies the original, unsquared resonance equation.
    distance = MISS_DISTANCE
    for x_root in x:
        if not math.isfinite(x_root):
            continue
        d_candidate = r0 * x_root - r * mu
        if d_candidate <= minimum_distance + distance_tolerance:
            continue
        if d_candidate > maximum_distance + distance_tolerance:
            continue

        new_r = math.sqrt(
            r * r + d_candidate * d_candidate + 2.0 * r * d_candidate * mu
        )
        if (
            new_r < r_inner - distance_tolerance
            or new_r > r_outer + distance_tolerance
        ):
            continue

        new_mu = (r * mu + d_candidate) / new_r
        new_v = v_inner + m * (new_r - r_inner)
        comov_nu = nu_rest * (
            1.0 - new_v / C_SPEED_OF_LIGHT * new_mu
        )
        if (
            not math.isfinite(comov_nu)
            or abs(comov_nu - nu_line) > 1.0e-7 * nu_line
        ):
            continue

        distance = min(distance, d_candidate)

    return distance


@njit(**njit_dict_no_parallel)
def calculate_projected_gradient_zero_distances(
    rpacket: RPacket,
    geometry: NumbaNonhomologousRadial1DGeometry,
    distance_boundary: float,
) -> tuple[float, float, int]:
    """Calculate forward projected-gradient zeros in the current shell.

    Returns up to two distances where the derivative of the projected fluid
    velocity along the packet trajectory changes sign. Missing distances are
    returned as ``MISS_DISTANCE``.
    """
    shell_id = rpacket.current_shell_id
    r = rpacket.r
    mu = rpacket.mu
    m = geometry.velocity_gradient[shell_id]
    q = geometry.v_outer[shell_id] - m * geometry.r_outer[shell_id]
    impact_parameter_squared = r * r * (1.0 - mu * mu)

    first_distance = MISS_DISTANCE
    second_distance = MISS_DISTANCE
    zero_count = 0

    if m == 0.0 or q == 0.0 or impact_parameter_squared == 0.0:
        return first_distance, second_distance, zero_count

    # For v(r) = m*r + q, the projected gradient along the ray is
    # m + q*b**2/r**3. Its zeros therefore satisfy r**3 = -q*b**2/m.
    turning_radius_cubed = -q * impact_parameter_squared / m
    if turning_radius_cubed <= 0.0:
        return first_distance, second_distance, zero_count

    turning_radius = turning_radius_cubed ** (1.0 / 3.0)
    shell_width = geometry.r_outer[shell_id] - geometry.r_inner[shell_id]
    distance_tolerance = CLOSE_LINE_THRESHOLD * shell_width
    impact_parameter = math.sqrt(impact_parameter_squared)

    # A zero at the impact parameter only touches zero and does not reverse
    # the monotonic line-list traversal direction.
    if turning_radius <= impact_parameter + distance_tolerance:
        return first_distance, second_distance, zero_count

    turning_x_squared = (
        turning_radius * turning_radius - impact_parameter_squared
    )
    turning_x = math.sqrt(turning_x_squared)
    initial_x = r * mu
    negative_x_distance = -turning_x - initial_x
    positive_x_distance = turning_x - initial_x

    if (
        negative_x_distance > distance_tolerance
        and negative_x_distance < distance_boundary - distance_tolerance
    ):
        first_distance = negative_x_distance
        zero_count = 1

    if (
        positive_x_distance > distance_tolerance
        and positive_x_distance < distance_boundary - distance_tolerance
    ):
        if zero_count == 0:
            first_distance = positive_x_distance
        else:
            second_distance = positive_x_distance
        zero_count += 1

    return first_distance, second_distance, zero_count


@njit(**njit_dict_no_parallel)
def calculate_comoving_frequency_nonhomologous(
    rpacket: RPacket,
    geometry: NumbaNonhomologousRadial1DGeometry,
    distance: float,
) -> float:
    """Calculate packet comoving frequency after a trajectory distance."""
    new_r = math.sqrt(
        rpacket.r * rpacket.r
        + distance * distance
        + 2.0 * rpacket.r * distance * rpacket.mu
    )
    new_mu = (rpacket.r * rpacket.mu + distance) / new_r
    new_v = geometry.get_velocity(new_r, rpacket.current_shell_id)
    return rpacket.nu * (1.0 - new_v / C_SPEED_OF_LIGHT * new_mu)


@njit(**njit_dict_no_parallel)
def get_line_id_range_nonhomologous(
    line_list_nu: npt.NDArray[np.float64],
    comov_nu_start: float,
    comov_nu_end: float,
) -> tuple[int, int, int]:
    """Return directional line-list bounds for one monotonic path interval.

    The transport line list is ordered by descending frequency. The returned
    range conservatively includes both interval endpoint frequencies; the
    distance bounds reject a resonance at the packet's interval start.
    """
    line_count = len(line_list_nu)
    frequency_tolerance = CLOSE_LINE_THRESHOLD * max(
        abs(comov_nu_start), abs(comov_nu_end)
    )
    minimum_frequency = (
        min(comov_nu_start, comov_nu_end) - frequency_tolerance
    )
    maximum_frequency = (
        max(comov_nu_start, comov_nu_end) + frequency_tolerance
    )

    # Find the first line with frequency <= maximum_frequency.
    lower_idx = 0
    upper_idx = line_count
    while lower_idx < upper_idx:
        middle_idx = (lower_idx + upper_idx) // 2
        if line_list_nu[middle_idx] > maximum_frequency:
            lower_idx = middle_idx + 1
        else:
            upper_idx = middle_idx
    maximum_frequency_line_id = lower_idx

    # Find the first line with frequency < minimum_frequency.
    lower_idx = 0
    upper_idx = line_count
    while lower_idx < upper_idx:
        middle_idx = (lower_idx + upper_idx) // 2
        if line_list_nu[middle_idx] >= minimum_frequency:
            lower_idx = middle_idx + 1
        else:
            upper_idx = middle_idx
    minimum_frequency_stop_line_id = lower_idx

    if comov_nu_end < comov_nu_start:
        return (
            maximum_frequency_line_id,
            minimum_frequency_stop_line_id,
            1,
        )

    return (
        minimum_frequency_stop_line_id - 1,
        maximum_frequency_line_id - 1,
        -1,
    )


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
