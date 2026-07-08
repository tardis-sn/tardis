"""Shared Monte Carlo transport loop utilities."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numba import njit, objmode, prange
from numba.np.ufunc.parallel import get_num_threads, get_thread_id
from numba.typed import List as TypedList

from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    create_estimators_bulk_list,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    create_estimators_line_list,
    init_estimators_line,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
    consolidate_vpacket_tracker,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.progress_bars import update_packets_pbar

if TYPE_CHECKING:
    from tardis.transport.montecarlo.configuration.base import (
        MonteCarloConfiguration,
    )
    from tardis.transport.montecarlo.packets.packet_collections import (
        PacketCollection,
    )


@njit(**njit_dict_no_parallel)
def make_r_packet(packet_collection: PacketCollection, packet_index: int):
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


@njit(**njit_dict_no_parallel)
def create_vpacket_collections(
    no_of_packets: int,
    spectrum_frequency_grid: np.ndarray,
    montecarlo_configuration: MonteCarloConfiguration,
    number_of_vpackets: int,
):
    """
    Create per-packet virtual packet collections.

    Parameters
    ----------
    no_of_packets : int
        Number of real packets.
    spectrum_frequency_grid : numpy.ndarray
        Frequency grid for virtual packet spectra.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.
    number_of_vpackets : int
        Number of virtual packets spawned per real packet interaction.

    Returns
    -------
    numba.typed.List
        Per-packet virtual packet collections.
    """
    vpacket_collections = TypedList()
    for i in range(no_of_packets):
        vpacket_collections.append(
            VPacketCollection(
                i,
                spectrum_frequency_grid,
                montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
                montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
                number_of_vpackets,
                montecarlo_configuration.TEMPORARY_V_PACKET_BINS,
            )
        )
    return vpacket_collections


@njit(**njit_dict_no_parallel)
def add_vpacket_collection_to_histogram(
    v_packets_energy_hist: np.ndarray,
    vpacket_collection: VPacketCollection,
    spectrum_frequency_grid: np.ndarray,
    delta_nu: float,
) -> None:
    """
    Add virtual packet collection energies to a spectrum histogram.

    Parameters
    ----------
    v_packets_energy_hist : numpy.ndarray
        Energy histogram to update in place.
    vpacket_collection : VPacketCollection
        Finalized virtual packet collection.
    spectrum_frequency_grid : numpy.ndarray
        Frequency grid for the histogram.
    delta_nu : float
        Frequency bin width.
    """
    v_packets_idx = np.floor(
        (vpacket_collection.nus - spectrum_frequency_grid[0]) / delta_nu
    ).astype(np.int64)

    for j, idx in enumerate(v_packets_idx):
        if (vpacket_collection.nus[j] < spectrum_frequency_grid[0]) or (
            vpacket_collection.nus[j] > spectrum_frequency_grid[-1]
        ):
            continue
        v_packets_energy_hist[idx] += vpacket_collection.energies[j]


@njit(**njit_dict_no_parallel)
def get_vpacket_tracker(
    vpacket_collections,
    spectrum_frequency_grid: np.ndarray,
    montecarlo_configuration: MonteCarloConfiguration,
) -> VPacketCollection:
    """
    Consolidate virtual packet tracking output when configured.

    Parameters
    ----------
    vpacket_collections : numba.typed.List
        Per-packet virtual packet collections.
    spectrum_frequency_grid : numpy.ndarray
        Frequency grid for virtual packet spectra.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.

    Returns
    -------
    VPacketCollection
        Consolidated tracker or an empty placeholder collection.
    """
    if montecarlo_configuration.ENABLE_VPACKET_TRACKING:
        return consolidate_vpacket_tracker(
            vpacket_collections,
            spectrum_frequency_grid,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
        )
    return VPacketCollection(
        -1,
        spectrum_frequency_grid,
        montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
        montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
        -1,
        1,
    )


@njit(**njit_dict)
def montecarlo_transport_with_vpackets(
    packet_collection: PacketCollection,
    geometry_state_numba,
    time_explosion: float,
    opacity_state_numba,
    montecarlo_configuration: MonteCarloConfiguration,
    spectrum_frequency_grid: np.ndarray,
    trackers,
    number_of_vpackets: int,
    show_progress_bars: bool,
    packet_propagation_function,
) -> tuple[
    np.ndarray,
    VPacketCollection,
    type,
    type,
]:
    """
    Run line-only Monte Carlo transport with virtual packet tracking.

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
    spectrum_frequency_grid : numpy.ndarray
        Frequency grid for virtual packet spectra.
    trackers
        Per-packet trackers.
    number_of_vpackets : int
        Number of virtual packets spawned per real packet interaction.
    show_progress_bars : bool
        Whether packet progress bars are enabled.
    packet_propagation_function
        Mode-specific packet propagation function.

    Returns
    -------
    tuple
        Virtual packet histogram, consolidated virtual packet tracker, bulk
        estimators, and line estimators.
    """
    no_of_packets = len(packet_collection.initial_nus)

    v_packets_energy_hist = np.zeros_like(spectrum_frequency_grid)
    delta_nu = spectrum_frequency_grid[1] - spectrum_frequency_grid[0]

    vpacket_collections = create_vpacket_collections(
        no_of_packets,
        spectrum_frequency_grid,
        montecarlo_configuration,
        number_of_vpackets,
    )

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
        vpacket_collection = vpacket_collections[packet_index]
        tracker = trackers[packet_index]

        packet_propagation_function(
            r_packet,
            geometry_state_numba,
            time_explosion,
            opacity_state_numba,
            estimators_bulk_thread,
            estimators_line_thread,
            vpacket_collection,
            tracker,
            montecarlo_configuration,
        )
        set_packet_collection_output(packet_collection, r_packet, i)

        tracker.finalize()
        vpacket_collection.finalize_arrays()

        add_vpacket_collection_to_histogram(
            v_packets_energy_hist,
            vpacket_collection,
            spectrum_frequency_grid,
            delta_nu,
        )

    for estimator_thread in estimators_bulk_list_thread:
        estimators_bulk.increment(estimator_thread)

    for estimator_thread in estimators_line_list_thread:
        estimators_line.increment(estimator_thread)

    vpacket_tracker = get_vpacket_tracker(
        vpacket_collections,
        spectrum_frequency_grid,
        montecarlo_configuration,
    )

    return (
        v_packets_energy_hist,
        vpacket_tracker,
        estimators_bulk,
        estimators_line,
    )
