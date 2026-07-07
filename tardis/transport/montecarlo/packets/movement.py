"""Shared radiative packet movement helpers."""

import numpy as np
from numba import njit

from tardis.transport.frame_transformations import (
    get_doppler_factor_from_velocity,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bulk,
)
from tardis.transport.montecarlo.packets.frame_transformations import (
    transform_packet_lab_to_comoving_frame,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)


@njit(**njit_dict_no_parallel)
def move_r_packet_with_velocity(
    r_packet: RPacket,
    distance: float,
    velocity: float,
    estimators_bulk: EstimatorsBulk,
    enable_full_relativity: bool,
) -> None:
    """
    Move a packet using a precomputed local ejecta velocity.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet object.
    distance : float
        Distance to move in cm.
    velocity : float
        Local ejecta velocity at the packet position before movement.
    estimators_bulk : EstimatorsBulk
        Cell-level bulk radiation field estimators.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.
    """
    doppler_factor = get_doppler_factor_from_velocity(
        velocity,
        r_packet.mu,
        enable_full_relativity,
    )

    if distance > 0.0:
        r = r_packet.r
        new_r = np.sqrt(
            r * r + distance * distance + 2.0 * r * distance * r_packet.mu
        )
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

        comov_nu = r_packet.nu * doppler_factor
        comov_energy = r_packet.energy * doppler_factor

        if enable_full_relativity:
            distance *= doppler_factor

        update_estimators_bulk(
            r_packet, distance, estimators_bulk, comov_nu, comov_energy
        )


@njit(**njit_dict_no_parallel)
def move_r_packet(
    r_packet: RPacket,
    distance: float,
    geometry,
    estimators_bulk: EstimatorsBulk,
    enable_full_relativity: bool,
) -> None:
    """
    Move a packet using local velocity from transport geometry.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet object.
    distance : float
        Distance to move in cm.
    geometry
        Geometry object exposing ``get_velocity(r, shell_id)``.
    estimators_bulk : EstimatorsBulk
        Cell-level bulk radiation field estimators.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.
    """
    comov_nu, comov_energy, _ = transform_packet_lab_to_comoving_frame(
        r_packet,
        geometry,
        enable_full_relativity,
    )

    if distance > 0.0:
        r = r_packet.r
        new_r = np.sqrt(
            r * r + distance * distance + 2.0 * r * distance * r_packet.mu
        )
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

        if enable_full_relativity:
            distance *= comov_nu / r_packet.nu

        update_estimators_bulk(
            r_packet, distance, estimators_bulk, comov_nu, comov_energy
        )


@njit(**njit_dict_no_parallel)
def increment_packet_cell_index(
    packet: RPacket, delta_shell: int, no_of_shells: int
) -> None:
    """
    Move a packet across a shell boundary and update its transport status.

    Parameters
    ----------
    packet : RPacket
        Radiative packet object.
    delta_shell : int
        Change in shell index, positive outward and negative inward.
    no_of_shells : int
        Number of shells in the transport grid.
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id
