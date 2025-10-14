import numpy as np
from numba import njit, objmode, prange
from numba.np.ufunc.parallel import get_num_threads, get_thread_id
from numba.typed import List

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import njit_dict
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
    VPacketCollection,
    consolidate_vpacket_tracker,
    initialize_last_interaction_tracker,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.progress_bars import update_packets_pbar
from tardis.transport.montecarlo.single_packet_loop import (
    single_packet_loop,
)


@njit(**njit_dict)
def montecarlo_main_loop(
    packet_collection: PacketCollection,
    geometry_state_numba: NumbaRadial1DGeometry,
    time_explosion: float,
    opacity_state_numba: OpacityStateNumba,
    montecarlo_configuration: MonteCarloConfiguration,
    estimators: RadiationFieldMCEstimators,
    spectrum_frequency_grid: np.ndarray,
    trackers: List,
    number_of_vpackets: int,
    show_progress_bars: bool,
):
    """
    Main loop of the Monte Carlo radiative transfer routine.

    This function generates the packet objects from the packet collection and propagates them through the ejecta,
    performing interactions and collecting statistics for the radiative
    transfer simulation.

    Parameters
    ----------
    packet_collection : PacketCollection
        Collection containing initial packet properties (positions, directions,
        frequencies, energies, and seeds)
    geometry_state : NumbaRadial1DGeometry
        Numba-compiled simulation geometry containing shell boundaries
        and velocity information
    time_explosion : float
        Time since explosion in seconds, used for relativistic calculations
    opacity_state : OpacityStateNumba
        Numba-compiled opacity state containing line opacities, continuum
        opacities, and atomic data required for interactions
    montecarlo_configuration : MonteCarloConfiguration
        Configuration object containing Monte Carlo simulation parameters
        and flags for various physics modules
    estimators : RadfieldMCEstimators
        Estimator object for collecting radiation field statistics
        during packet propagation
    spectrum_frequency_grid : numpy.ndarray
        Frequency grid array for virtual packet spectrum calculation
    rpacket_trackers : numba.typed.List
        List of packet trackers for detailed packet interaction logging
    number_of_vpackets : int
        Number of virtual packets to spawn per real packet interaction
    show_progress_bars : bool
        Flag to enable/disable progress bar updates during simulation

    Returns
    -------
    tuple
        A tuple containing:
        - v_packets_energy_hist : numpy.ndarray
            Energy histogram of virtual packets binned by frequency
        - last_interaction_tracker : LastInteractionTracker
            Object tracking the last interaction properties for each packet
        - vpacket_tracker : VPacketCollection
            Consolidated virtual packet collection containing all virtual
            packet information from the simulation
    """
    no_of_packets = len(packet_collection.initial_nus)

    v_packets_energy_hist = np.zeros_like(spectrum_frequency_grid)
    delta_nu = spectrum_frequency_grid[1] - spectrum_frequency_grid[0]

    # Pre-allocate a list of vpacket collections for later storage
    vpacket_collections = List()
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

    # Get the ID of the main thread and the number of threads
    main_thread_id = get_thread_id()
    n_threads = get_num_threads()

    # betting get thread_id goes from 0 to num threads
    # Note that get_thread_id() returns values from 0 to n_threads-1,
    # so we iterate from 0 to n_threads-1 to create the estimator_list
    estimator_list = estimators.create_estimator_list(n_threads)

    for i in prange(no_of_packets):
        thread_id = get_thread_id()
        if show_progress_bars:
            if thread_id == main_thread_id:
                with objmode:
                    update_amount = 1 * n_threads
                    update_packets_pbar(
                        update_amount,
                        no_of_packets,
                    )

        r_packet = RPacket(
            packet_collection.initial_radii[i],
            packet_collection.initial_mus[i],
            packet_collection.initial_nus[i],
            packet_collection.initial_energies[i],
            packet_collection.packet_seeds[i],
            i,
        )
        # Seed the random number generator
        np.random.seed(r_packet.seed)

        # Get the local estimators for this thread
        local_estimators = estimator_list[thread_id]

        # Get the local v_packet_collection for this thread
        vpacket_collection = vpacket_collections[i]
        # RPacket Tracker for this thread
        tracker = trackers[i]

        loop = single_packet_loop(
            r_packet,
            geometry_state_numba,
            time_explosion,
            opacity_state_numba,
            local_estimators,
            vpacket_collection,
            tracker,
            montecarlo_configuration,
        )
        packet_collection.output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            packet_collection.output_energies[i] = -r_packet.energy

        elif r_packet.status == PacketStatus.EMITTED:
            packet_collection.output_energies[i] = r_packet.energy

        # Finalize the tracker (e.g. trim arrays to actual size)
        tracker.finalize()
        
        # Finalize the vpacket collection to trim arrays to actual size
        vpacket_collection.finalize_arrays()

        v_packets_idx = np.floor(
            (vpacket_collection.nus - spectrum_frequency_grid[0]) / delta_nu
        ).astype(np.int64)

        for j, idx in enumerate(v_packets_idx):
            if (vpacket_collection.nus[j] < spectrum_frequency_grid[0]) or (
                vpacket_collection.nus[j] > spectrum_frequency_grid[-1]
            ):
                continue
            v_packets_energy_hist[idx] += vpacket_collection.energies[j]

    for sub_estimator in estimator_list:
        estimators.increment(sub_estimator)

    if montecarlo_configuration.ENABLE_VPACKET_TRACKING:
        vpacket_tracker = consolidate_vpacket_tracker(
            vpacket_collections,
            spectrum_frequency_grid,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
        )
    else:
        vpacket_tracker = VPacketCollection(
            -1,
            spectrum_frequency_grid,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
            -1,
            1,
        )

    return (
        v_packets_energy_hist,
        vpacket_tracker,
    )
