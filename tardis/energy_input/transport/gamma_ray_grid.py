import numpy as np
from numba import njit
from numpy.typing import NDArray

from tardis.energy_input.transport.GXPacket import GXPacket
from tardis.energy_input.util import (
    C_CGS,
    doppler_factor_3d,
    solve_quadratic_equation,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


@njit(**njit_dict_no_parallel)
def calculate_distance_radial(
    photon: GXPacket, r_inner: float, r_outer: float
) -> tuple[float, int]:
    """Calculate the 3D distance to a gamma-ray shell boundary.

    Parameters
    ----------
    photon : GXPhoton object
    r_inner : float
    r_outer : float

    Returns
    -------
    distance : float
        Distance to the selected shell boundary in centimeters.
    shell_change : int
        Change in radial-shell index at the selected boundary.

    Notes
    -----
    A zero inner radius does not define a boundary. Packets in that central
    shell continue through the origin to the far-side outer boundary.
    """
    # solve the quadratic distance equation for the inner and
    # outer shell boundaries
    inner_1 = -1.0
    inner_2 = -1.0
    if r_inner > 0.0:
        inner_1, inner_2 = solve_quadratic_equation(
            photon.location, photon.direction, r_inner
        )
    outer_1, outer_2 = solve_quadratic_equation(
        photon.location, photon.direction, r_outer
    )

    final_position_outer_1 = np.ascontiguousarray(
        photon.location + photon.direction * outer_1
    )
    final_position_outer_2 = np.ascontiguousarray(
        photon.location + photon.direction * outer_2
    )

    # Ensure photon.direction is contiguous for dot product operations
    direction_contiguous = np.ascontiguousarray(photon.direction)

    if r_inner > 0.0:
        final_position_inner_1 = np.ascontiguousarray(
            photon.location + photon.direction * inner_1
        )
        final_position_inner_2 = np.ascontiguousarray(
            photon.location + photon.direction * inner_2
        )
        if np.dot(final_position_inner_1, direction_contiguous) > 0:
            inner_1 = -1
        if np.dot(final_position_inner_2, direction_contiguous) > 0:
            inner_2 = -1
    if np.dot(final_position_outer_1, direction_contiguous) < 0:
        outer_1 = -1
    if np.dot(final_position_outer_2, direction_contiguous) < 0:
        outer_2 = -1

    distances = np.array([inner_1, inner_2, outer_1, outer_2])

    # the correct distance is the shortest positive distance
    distance_list = [i for i in distances if i > 0]

    if not distance_list:
        print(photon.get_location_r() - r_inner)
        print(photon.get_location_r() - r_outer)
        print(photon.get_location_r())
        print(photon.location, photon.direction, r_inner, r_outer)
        print(distances)
        print(photon.shell)
        raise ValueError("No root found for distance calculation!")

    shortest = min(distance_list)
    shell_change = 1

    if shortest in (inner_1, inner_2):
        shell_change = -1

    return shortest, shell_change


@njit(**njit_dict_no_parallel)
def distance_trace(
    photon: GXPacket,
    inner_velocity: NDArray[np.float64],
    outer_velocity: NDArray[np.float64],
    total_opacity: float,
    current_time: float,
    next_time: float,
) -> tuple[float, float, float, int]:
    """Trace distances to the next gamma-ray event.

    Calculate the distances to an interaction, shell boundary, and time-bin
    boundary.

    Parameters
    ----------
    photon : GXPhoton object
    inner_velocity : One dimensional Numpy array, dtype float
    outer_velocity : One dimensional Numpy array, dtype float
    total_opacity : float
    current_time : float
    next_time : float

    Returns
    -------
    distance_interaction : float
    distance_boundary : float
    distance_time : float
    shell_change : int

    Notes
    -----
    The first high-energy transport shell extends to the origin. Its geometric
    inner boundary is ignored so that inward packets cross the center without
    leaving the radial grid.
    """
    r_inner = 0.0
    if photon.shell > 0:
        r_inner = inner_velocity[photon.shell] * current_time

    distance_boundary, shell_change = calculate_distance_radial(
        photon,
        r_inner,
        outer_velocity[photon.shell] * current_time,
    )

    distance_interaction = photon.tau / total_opacity
    distance_time = (next_time - photon.time_start) * C_CGS
    return distance_interaction, distance_boundary, distance_time, shell_change


@njit(**njit_dict_no_parallel)
def move_packet(packet: GXPacket, distance: float) -> GXPacket:
    """Move a packet a distance along its direction vector.

    Parameters
    ----------
    packet : GXPacket object
    distance : float

    Returns
    -------
    packet : GXPacket object

    """
    location_old = packet.location
    direction = packet.direction

    location_new = location_old + distance * direction

    packet.location = location_new

    doppler_factor = doppler_factor_3d(
        packet.direction, packet.location, packet.time_start
    )

    packet.nu_cmf = packet.nu_rf * doppler_factor
    packet.energy_cmf = packet.energy_rf * doppler_factor

    return packet
