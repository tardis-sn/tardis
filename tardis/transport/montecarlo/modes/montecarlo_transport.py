"""Shared Monte Carlo transport loop utilities."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numba import njit, objmode, prange
from numba.np.ufunc.parallel import get_num_threads, get_thread_id
from numba.typed import List as TypedList

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    create_estimators_bulk_list,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
    create_estimators_line_list,
    init_estimators_line,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)
from tardis.transport.montecarlo.progress_bars import update_packets_pbar


@njit(**njit_dict_no_parallel)
def make_r_packet(
    packet_collection: PacketCollection, packet_index: int
) -> RPacket:
    """
    Build and seed a radiative packet from a packet collection.

    Parameters
    ----------
    packet_collection : PacketCollection
        Packet collection containing the initial packet arrays.
    packet_index : int
        Packet index.

    Returns
    -------
    RPacket
        Radiative packet initialized from the packet collection.
    """
    r_packet = RPacket(
        packet_collection.initial_radii[packet_index],
        packet_collection.initial_mus[packet_index],
        packet_collection.initial_nus[packet_index],
        packet_collection.initial_energies[packet_index],
        packet_collection.packet_seeds[packet_index],
        packet_index,
    )
    np.random.seed(r_packet.seed)
    return r_packet


@njit(**njit_dict_no_parallel)
def set_packet_collection_output(
    packet_collection: PacketCollection, r_packet: RPacket, packet_index: int
) -> None:
    """
    Store final packet frequency and signed energy in a packet collection.

    Parameters
    ----------
    packet_collection : PacketCollection
        Packet collection to update.
    r_packet : RPacket
        Packet after propagation.
    packet_index : int
        Packet index.
    """
    packet_collection.output_nus[packet_index] = r_packet.nu

    if r_packet.status == PacketStatus.REABSORBED:
        packet_collection.output_energies[packet_index] = -r_packet.energy
    elif r_packet.status == PacketStatus.EMITTED:
        packet_collection.output_energies[packet_index] = r_packet.energy


@njit(**njit_dict_no_parallel)
def update_packet_progress(
    show_progress_bars: bool,
    thread_id: int,
    main_thread_id: int,
    n_threads: int,
    no_of_packets: int,
) -> None:
    """
    Update the packet progress bar from the main Numba worker thread.

    Parameters
    ----------
    show_progress_bars : bool
        Whether progress bars are enabled.
    thread_id : int
        Current worker thread id.
    main_thread_id : int
        Thread id selected to perform progress updates.
    n_threads : int
        Number of Numba worker threads.
    no_of_packets : int
        Total number of packets.
    """
    if show_progress_bars:
        if thread_id == main_thread_id:
            with objmode:
                update_packets_pbar(n_threads, no_of_packets)


@njit(**njit_dict)
def montecarlo_transport(
    packet_collection: PacketCollection,
    geometry_state_numba: (
        NumbaHomologousRadial1DGeometry | NumbaRadial1DGeometry
    ),
    time_explosion: float,
    opacity_state_numba: OpacityStateNumba,
    montecarlo_configuration: MonteCarloConfiguration,
    trackers: (TypedList[TrackerFull] | TypedList[TrackerLastInteraction]),
    show_progress_bars: bool,
    packet_propagation_function: Callable[..., None],
) -> tuple[EstimatorsBulk, EstimatorsLine]:
    """
    Run real-packet Monte Carlo transport and record packet histories.

    Parameters
    ----------
    packet_collection : PacketCollection
        Packet collection containing packet input and output arrays.
    geometry_state_numba
        Numba geometry object for the transport mode.
    time_explosion : float
        Time since explosion in seconds. Non-homologous mode accepts but does
        not use this value.
    opacity_state_numba
        Numba opacity state.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.
    trackers
        Per-packet trackers.
    show_progress_bars : bool
        Whether packet progress bars are enabled.
    packet_propagation_function
        Mode-specific packet propagation function.

    Returns
    -------
    tuple
        Bulk and line estimators.
    """
    no_of_packets = len(packet_collection.initial_nus)

    main_thread_id = get_thread_id()
    n_threads = get_num_threads()

    n_lines_by_n_cells_tuple = opacity_state_numba.tau_sobolev.shape
    n_cells = len(geometry_state_numba.r_inner)
    estimators_bulk = init_estimators_bulk(n_cells)
    estimators_line = init_estimators_line(n_lines_by_n_cells_tuple)

    estimators_bulk_list_thread = create_estimators_bulk_list(
        n_cells, n_threads
    )
    estimators_line_list_thread = create_estimators_line_list(
        n_lines_by_n_cells_tuple, n_threads
    )

    for i in prange(no_of_packets):
        packet_index = np.int64(i)
        thread_id = get_thread_id()
        update_packet_progress(
            show_progress_bars,
            thread_id,
            main_thread_id,
            n_threads,
            no_of_packets,
        )

        r_packet = make_r_packet(packet_collection, packet_index)
        estimators_bulk_thread = estimators_bulk_list_thread[thread_id]
        estimators_line_thread = estimators_line_list_thread[thread_id]
        tracker = trackers[packet_index]

        packet_propagation_function(
            r_packet,
            geometry_state_numba,
            time_explosion,
            opacity_state_numba,
            estimators_bulk_thread,
            estimators_line_thread,
            tracker,
            montecarlo_configuration,
        )
        set_packet_collection_output(packet_collection, r_packet, i)

        tracker.finalize()

    for estimator_thread in estimators_bulk_list_thread:
        estimators_bulk.increment(estimator_thread)

    for estimator_thread in estimators_line_list_thread:
        estimators_line.increment(estimator_thread)

    return estimators_bulk, estimators_line
