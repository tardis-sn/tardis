import math

import numpy as np
import numpy.typing as npt
from numba import njit

from tardis.model.geometry.radial1d import (
    NumbaRadial1DGeometry,
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
    solve_resonance_quartic,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket
from tardis.transport.montecarlo.utils import MonteCarloException

RESONANCE_FREQUENCY_RELATIVE_TOLERANCE = 1.0e-7


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
    r_packet: RPacket,
    comov_nu: float,
    is_last_line: bool,
    nu_line: float,
    time_explosion: float,
    enable_full_relativity: bool,
) -> float:
    """
    Calculate the homologous-flow distance to the next line resonance.

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Radiative packet being propagated.
    comov_nu : float
        Comoving frequency at the current packet position.
    is_last_line : bool
        Whether the packet has reached the end of the line list.
    nu_line : float
        Rest frequency of the target line.
    time_explosion : float
        Time since explosion in seconds.
    enable_full_relativity : bool
        Whether to use the full-relativity distance calculation.

    Returns
    -------
    float
        Distance to the line resonance in centimeters.
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
def calculate_packet_velocity_properties(
    rpacket: RPacket,
    geometry: NumbaRadial1DGeometry,
    distance: float,
) -> tuple[float, float, float, float]:
    """Calculate packet position and local velocity-field properties.

    Parameters
    ----------
    rpacket : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Radiative packet at the start of its trajectory segment.
    geometry : tardis.model.geometry.radial1d.NumbaRadial1DGeometry
        Radial geometry containing the piecewise-linear velocity field.
    distance : float
        Distance traveled along the packet trajectory in centimeters.

    Returns
    -------
    tuple[float, float, float, float]
        Radius, direction cosine, radial velocity, and projected velocity
        gradient after traveling ``distance``.
    """
    radius = math.sqrt(
        rpacket.r * rpacket.r
        + distance * distance
        + 2.0 * rpacket.r * distance * rpacket.mu
    )
    direction_cosine = (rpacket.r * rpacket.mu + distance) / radius
    velocity = geometry.get_velocity(radius, rpacket.current_shell_id)
    radial_velocity_gradient = geometry.velocity_gradient[
        rpacket.current_shell_id
    ]
    projected_velocity_gradient = (
        direction_cosine**2 * radial_velocity_gradient
        + (1.0 - direction_cosine**2) * velocity / radius
    )
    return radius, direction_cosine, velocity, projected_velocity_gradient


@njit(fastmath=False, error_model="numpy", parallel=False)
def calculate_distance_line_nonhomologous(
    rpacket: RPacket,
    geometry: NumbaRadial1DGeometry,
    nu_line: float,
    minimum_distance: float = 0.0,
    maximum_distance: float = MISS_DISTANCE,
) -> float:
    """
    Calculate the distance to a line resonance in a radial velocity field.

    Candidate roots must lie after ``minimum_distance`` and no farther than
    ``maximum_distance`` along the current packet trajectory.

    Parameters
    ----------
    rpacket : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Radiative packet being propagated.
    geometry : tardis.model.geometry.radial1d.NumbaRadial1DGeometry
        Radial geometry containing the piecewise-linear velocity field.
    nu_line : float
        Rest frequency of the target line.
    minimum_distance : float, optional
        Exclusive lower bound for candidate resonance distances.
    maximum_distance : float, optional
        Inclusive upper bound for candidate resonance distances.

    Returns
    -------
    float
        Distance to the line resonance in centimeters.
    """
    r_inner = geometry.r_inner[rpacket.current_shell_id]
    r_outer = geometry.r_outer[rpacket.current_shell_id]
    v_inner = geometry.v_inner[rpacket.current_shell_id]
    v_outer = geometry.v_outer[rpacket.current_shell_id]

    radius = rpacket.r
    velocity_gradient = geometry.velocity_gradient[rpacket.current_shell_id]
    rest_frequency = rpacket.nu
    direction_cosine = rpacket.mu

    target_projected_velocity = C_SPEED_OF_LIGHT * (
        1.0 - nu_line / rest_frequency
    )
    transverse_direction_fraction = 1.0 - direction_cosine**2
    velocity_intercept = v_outer - velocity_gradient * r_outer

    # Characteristic scales for non-dimensionalization
    shell_width = r_outer - r_inner
    shell_velocity_difference = v_outer - v_inner
    distance_tolerance = CLOSE_LINE_THRESHOLD * shell_width
    if shell_velocity_difference == 0.0:
        impact_parameter_squared = (
            radius * radius * transverse_direction_fraction
        )
        if (
            velocity_intercept != 0.0
            and target_projected_velocity**2 < velocity_intercept**2
        ):
            projected_position_squared = (
                target_projected_velocity**2
                * impact_parameter_squared
                / (velocity_intercept**2 - target_projected_velocity**2)
            )
            scaled_projected_position = (
                math.copysign(
                    math.sqrt(projected_position_squared),
                    target_projected_velocity / velocity_intercept,
                )
                / shell_width
            )
            projected_position_roots = (
                scaled_projected_position,
                math.nan,
                math.nan,
                math.nan,
            )
        else:
            projected_position_roots = (
                math.nan,
                math.nan,
                math.nan,
                math.nan,
            )
    else:
        # Dimensionless quantities improve quartic-solver floating-point accuracy.
        scaled_radius = radius / shell_width
        scaled_target_projected_velocity = (
            target_projected_velocity / shell_velocity_difference
        )
        scaled_velocity_intercept = (
            velocity_intercept / shell_velocity_difference
        )
        scaled_impact_parameter_squared = (
            scaled_radius**2 * transverse_direction_fraction
        )

        # Obtain roots of the quartic polynomial for the dimensionless
        # x = (d_line + r_i \mu_i) / r0.
        projected_position_roots = solve_resonance_quartic(
            scaled_target_projected_velocity,
            scaled_impact_parameter_squared,
            scaled_velocity_intercept,
        )

    # Convert each dimensionless root to a distance: d = r0*x_i - r*mu.
    # Select the nearest root that satisfies the original, unsquared resonance equation.
    distance = MISS_DISTANCE
    for projected_position_root in projected_position_roots:
        if not math.isfinite(projected_position_root):
            continue
        candidate_distance = (
            shell_width * projected_position_root - radius * direction_cosine
        )
        if candidate_distance <= minimum_distance + distance_tolerance:
            continue
        if candidate_distance > maximum_distance + distance_tolerance:
            continue

        (
            resonance_radius,
            resonance_direction_cosine,
            velocity_at_resonance,
            _,
        ) = calculate_packet_velocity_properties(
            rpacket, geometry, candidate_distance
        )
        if (
            resonance_radius < r_inner - distance_tolerance
            or resonance_radius > r_outer + distance_tolerance
        ):
            continue

        comoving_frequency_at_resonance = rest_frequency * (
            1.0
            - velocity_at_resonance
            / C_SPEED_OF_LIGHT
            * resonance_direction_cosine
        )
        if (
            not math.isfinite(comoving_frequency_at_resonance)
            or abs(comoving_frequency_at_resonance - nu_line)
            > RESONANCE_FREQUENCY_RELATIVE_TOLERANCE * nu_line
        ):
            continue

        distance = min(distance, candidate_distance)

    return distance


@njit(**njit_dict_no_parallel)
def calculate_projected_gradient_zero_distances(
    rpacket: RPacket,
    geometry: NumbaRadial1DGeometry,
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
    geometry: NumbaRadial1DGeometry,
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
    minimum_frequency = min(comov_nu_start, comov_nu_end) - frequency_tolerance
    maximum_frequency = max(comov_nu_start, comov_nu_end) + frequency_tolerance

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
    nu_line: float,
    nu: float,
    time_explosion: float,
    r_packet: RPacket,
) -> float:
    """Calculate a fully relativistic homologous line-resonance distance."""
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
