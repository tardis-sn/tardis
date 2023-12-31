import numpy as np
from numba import njit, objmode, prange
from numba.np.ufunc.parallel import get_num_threads, get_thread_id
from numba.typed import List

from tardis.montecarlo import montecarlo_configuration
from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis.montecarlo.estimators.estimator_statistics import (
    EstimatorStatistics,
)
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    NumbaModel,
    RPacketTracker,
)
from tardis.montecarlo.montecarlo_numba.packet_collections import (
    VPacketCollection,
    consolidate_vpacket_tracker,
    initialize_last_interaction_tracker,
)
from tardis.montecarlo.montecarlo_numba.r_packet import (
    PacketStatus,
    RPacket,
)
from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop,
)
from tardis.util.base import update_packet_pbar


@njit(**njit_dict)
def montecarlo_main_loop(
    packet_collection,
    geometry_state,
    numba_model,
    opacity_state,
    estimators,
    spectrum_frequency,
    number_of_vpackets,
    iteration,
    show_progress_bars,
    total_iterations,
):
    """
    This is the main loop of the MonteCarlo routine that generates packets
    and sends them through the ejecta.

    Parameters
    ----------
    packet_collection : PacketCollection
    numba_radial_1d_geometry : NumbaRadial1DGeometry
    numba_model : NumbaModel
    opacity_state : OpacityState
    estimators : NumbaEstimators
    spectrum_frequency : astropy.units.Quantity
        frequency binspas
    number_of_vpackets : int
        VPackets released per interaction
    packet_seeds : numpy.array
    virtual_packet_logging : bool
        Option to enable virtual packet logging.
    """
    no_of_packets = len(packet_collection.initial_nus)

    last_interaction_tracker = initialize_last_interaction_tracker(
        no_of_packets
    )

    v_packets_energy_hist = np.zeros_like(spectrum_frequency)
    delta_nu = spectrum_frequency[1] - spectrum_frequency[0]

    # Pre-allocate a list of vpacket collections for later storage
    vpacket_collections = List()
    # Configuring the Tracking for R_Packets
    rpacket_trackers = List()
    for i in range(no_of_packets):
        vpacket_collections.append(
            VPacketCollection(
                i,
                spectrum_frequency,
                montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
                montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
                number_of_vpackets,
                montecarlo_configuration.TEMPORARY_V_PACKET_BINS,
            )
        )
        rpacket_trackers.append(RPacketTracker())

    # Get the ID of the main thread and the number of threads
    main_thread_id = get_thread_id()
    n_threads = get_num_threads()

    estimator_list = List()
    for i in range(
        n_threads
    ):  # betting get thread_id goes from 0 to num threads
        # Note that get_thread_id() returns values from 0 to n_threads-1,
        # so we iterate from 0 to n_threads-1 to create the estimator_list
        estimator_list.append(
            EstimatorStatistics(
                np.copy(estimators.j_estimator),
                np.copy(estimators.nu_bar_estimator),
                np.copy(estimators.j_blue_estimator),
                np.copy(estimators.Edotlu_estimator),
                np.copy(estimators.photo_ion_estimator),
                np.copy(estimators.stim_recomb_estimator),
                np.copy(estimators.bf_heating_estimator),
                np.copy(estimators.stim_recomb_cooling_estimator),
                np.copy(estimators.photo_ion_estimator_statistics),
            )
        )

    for i in prange(no_of_packets):
        thread_id = get_thread_id()
        if show_progress_bars:
            if thread_id == main_thread_id:
                with objmode:
                    update_amount = 1 * n_threads
                    update_packet_pbar(
                        update_amount,
                        current_iteration=iteration,
                        no_of_packets=no_of_packets,
                        total_iterations=total_iterations,
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
        rpacket_tracker = rpacket_trackers[i]

        loop = single_packet_loop(
            r_packet,
            geometry_state,
            numba_model,
            opacity_state,
            local_estimators,
            vpacket_collection,
            rpacket_tracker,
        )
        packet_collection.output_nus[i] = r_packet.nu

        last_interaction_tracker.update_last_interaction(r_packet, i)

        if r_packet.status == PacketStatus.REABSORBED:
            packet_collection.output_energies[i] = -r_packet.energy
            last_interaction_tracker.types[i] = r_packet.last_interaction_type
        elif r_packet.status == PacketStatus.EMITTED:
            packet_collection.output_energies[i] = r_packet.energy
            last_interaction_tracker.types[i] = r_packet.last_interaction_type

        vpacket_collection.finalize_arrays()

        v_packets_idx = np.floor(
            (vpacket_collection.nus - spectrum_frequency[0]) / delta_nu
        ).astype(np.int64)

        for j, idx in enumerate(v_packets_idx):
            if (vpacket_collection.nus[j] < spectrum_frequency[0]) or (
                vpacket_collection.nus[j] > spectrum_frequency[-1]
            ):
                continue
            v_packets_energy_hist[idx] += vpacket_collection.energies[j]

    for sub_estimator in estimator_list:
        estimators.increment(sub_estimator)

    if montecarlo_configuration.ENABLE_VPACKET_TRACKING:
        vpacket_tracker = consolidate_vpacket_tracker(
            vpacket_collections, spectrum_frequency
        )
    else:
        vpacket_tracker = VPacketCollection(
            -1,
            spectrum_frequency,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
            -1,
            1,
        )

    if montecarlo_configuration.ENABLE_RPACKET_TRACKING:
        for rpacket_tracker in rpacket_trackers:
            rpacket_tracker.finalize_array()

    return (
        v_packets_energy_hist,
        last_interaction_tracker,
        vpacket_tracker,
        rpacket_trackers,
    )
