import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.energy_input.util import (
    doppler_factor_3d,
    solve_quadratic_equation,
    C_CGS,
)


@njit(**njit_dict_no_parallel)
def calculate_distance_radial(photon, r_inner, r_outer):
    """
    Calculates 3D distance to shell from gamma ray position

    Parameters
    ----------
    photon : GXPhoton object
    r_inner : float
    r_outer : float

    Returns
    -------
    distance : float

    """

    # solve the quadratic distance equation for the inner and
    # outer shell boundaries
    inner_1, inner_2 = solve_quadratic_equation(
        photon.location, photon.direction, r_inner
    )
    outer_1, outer_2 = solve_quadratic_equation(
        photon.location, photon.direction, r_outer
    )

    final_position_inner_1 = photon.location + photon.direction * inner_1
    final_position_inner_2 = photon.location + photon.direction * inner_2
    final_position_outer_1 = photon.location + photon.direction * outer_1
    final_position_outer_2 = photon.location + photon.direction * outer_2

    if np.dot(final_position_inner_1, photon.direction) > 0:
        inner_1 = -1
    if np.dot(final_position_inner_2, photon.direction) > 0:
        inner_2 = -1
    if np.dot(final_position_outer_1, photon.direction) < 0:
        outer_1 = -1
    if np.dot(final_position_outer_2, photon.direction) < 0:
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

    if shortest == (inner_1 or inner_2):
        shell_change = -1

    return shortest, shell_change


@njit(**njit_dict_no_parallel)
def distance_trace(
    photon,
    inner_velocity,
    outer_velocity,
    total_opacity,
    current_time,
    next_time,
):
    """
    Traces distance traveled by gamma ray and finds distance to
    next interaction and boundary

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
    """
    distance_boundary, shell_change = calculate_distance_radial(
        photon,
        inner_velocity[photon.shell] * current_time,
        outer_velocity[photon.shell] * current_time,
    )

    distance_interaction = photon.tau / total_opacity
    distance_time = (next_time - photon.time_current) * C_CGS
    return distance_interaction, distance_boundary, distance_time, shell_change


@njit(**njit_dict_no_parallel)
def move_packet(packet, distance):
    """
    Moves packet a distance along its direction vector

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
        packet.direction, packet.location, packet.time_current
    )

    packet.nu_cmf = packet.nu_rf * doppler_factor
    packet.energy_cmf = packet.energy_rf * doppler_factor

    return packet
