"""Packet movement routines shared by Monte Carlo transport modes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from tardis.transport.frame_transformations import get_doppler_factor
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bulk,
)
from tardis.transport.montecarlo.packets.radiative_packet import PacketStatus

if TYPE_CHECKING:
    from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
    from tardis.model.geometry.radial1d_nonhomologous import (
        NumbaNonhomologousRadial1DGeometry,
    )
    from tardis.transport.montecarlo.estimators.estimators_bulk import (
        EstimatorsBulk,
    )
    from tardis.transport.montecarlo.packets.radiative_packet import RPacket


@njit(**njit_dict_no_parallel)
def move_r_packet(
    r_packet: RPacket,
    distance: float,
    geometry: NumbaRadial1DGeometry | NumbaNonhomologousRadial1DGeometry,
    estimators_bulk: EstimatorsBulk,
    enable_full_relativity: bool,
    full_relativity_supported: bool = True,
) -> None:
    """
    Move a radiative packet using velocity from the model geometry.

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Radiative packet to move.
    distance : float
        Lab-frame distance traveled by the packet [cm].
    geometry : NumbaRadial1DGeometry or NumbaNonhomologousRadial1DGeometry
        Geometry object that provides local packet-frame velocity.
    estimators_bulk : tardis.transport.montecarlo.estimators.estimators_bulk.EstimatorsBulk
        Cell-level bulk radiation field estimators to update in place.
    enable_full_relativity : bool
        Whether to apply full-relativity Doppler and path-length corrections.
    full_relativity_supported : bool, optional
        Whether the supplied geometry supports full relativity in this
        transport mode, by default True.
    """
    if enable_full_relativity and not full_relativity_supported:
        raise NotImplementedError(
            "Full relativity not implemented for non-homologous mode."
        )

    velocity = geometry.get_velocity(r_packet.r, r_packet.current_shell_id)
    doppler_factor = get_doppler_factor(
        velocity, r_packet.mu, enable_full_relativity
    )

    r = r_packet.r
    if distance > 0.0:
        new_r = np.sqrt(
            r * r + distance * distance + 2.0 * r * distance * r_packet.mu
        )
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

        comov_nu = r_packet.nu * doppler_factor
        comov_energy = r_packet.energy * doppler_factor

        # Account for length contraction
        if enable_full_relativity:
            distance *= doppler_factor

        update_estimators_bulk(
            r_packet, distance, estimators_bulk, comov_nu, comov_energy
        )


@njit(**njit_dict_no_parallel)
def move_packet_across_shell_boundary(
    packet: RPacket, delta_shell: int, no_of_shells: int
) -> None:
    """
    Move packet across a shell boundary and update packet status.

    Parameters
    ----------
    packet : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Packet object.
    delta_shell : int
        Change in shell index (+1 if moving outward or -1 if moving inward).
    no_of_shells : int
        Number of shells in the simulation.
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id
