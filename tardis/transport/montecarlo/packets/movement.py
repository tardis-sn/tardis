"""Packet movement routines shared by Monte Carlo transport modes."""

from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bulk,
)

if TYPE_CHECKING:
    from tardis.transport.montecarlo.estimators.estimators_bulk import (
        EstimatorsBulk,
    )
    from tardis.transport.montecarlo.packets.radiative_packet import RPacket


@njit(**njit_dict_no_parallel)
def move_r_packet_core(
    r_packet: RPacket,
    distance: float,
    doppler_factor: float,
    enable_full_relativity: bool,
    estimators_bulk: EstimatorsBulk,
) -> None:
    """
    Move a radiative packet and update bulk radiation field estimators.

    The packet radius and directional cosine are updated in place after
    transport through `distance` in the lab frame. Bulk estimators are updated
    with comoving-frame frequency and energy; when full relativity is enabled,
    the estimator path length is also transformed with the supplied Doppler
    factor.

    Parameters
    ----------
    r_packet : tardis.transport.montecarlo.packets.radiative_packet.RPacket
        Radiative packet to move.
    distance : float
        Lab-frame distance traveled by the packet [cm].
    doppler_factor : float
        Doppler factor used to transform packet frequency and energy to the
        comoving frame.
    enable_full_relativity : bool
        Whether to apply full-relativity path-length correction to estimator
        updates.
    estimators_bulk : tardis.transport.montecarlo.estimators.estimators_bulk.EstimatorsBulk
        Cell-level bulk radiation field estimators to update in place.
    """
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
