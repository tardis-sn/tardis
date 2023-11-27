from numba import prange, njit, objmode
from numba.np.ufunc.parallel import get_thread_id, get_num_threads

import numpy as np
from tardis.montecarlo.montecarlo_numba.estimators import Estimators
from tardis.montecarlo.montecarlo_numba.packet_collections import (
    VPacketCollection,
    initialize_last_interaction_tracker,
)


from tardis.montecarlo.montecarlo_numba.r_packet import (
    RPacket,
    PacketStatus,
)

from tardis.montecarlo.montecarlo_numba.numba_interface import (
    RPacketTracker,
    NumbaModel,
)

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop,
)
from tardis.montecarlo.montecarlo_numba import njit_dict
from numba.typed import List
from tardis.util.base import (
    update_iterations_pbar,
    update_packet_pbar,
    refresh_packet_pbar,
)


def montecarlo_radial1d(
    transport_state,
    time_explosion,
    iteration,
    total_iterations,
    show_progress_bars=False,
):
    packet_collection = transport_state.packet_collection
    estimators = transport_state.estimators
    geometry_state = transport_state.geometry_state
    opacity_state = transport_state.opacity_state
    numba_model = NumbaModel(
        time_explosion.to("s").value,
    )

    number_of_vpackets = montecarlo_configuration.number_of_vpackets

    (
        v_packets_energy_hist,
        last_interaction_tracker,
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_initial_mus,
        virt_packet_initial_rs,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
        virt_packet_last_line_interaction_shell_id,
        rpacket_trackers,
    ) = montecarlo_main_loop(
        packet_collection,
        geometry_state,
        numba_model,
        opacity_state,
        estimators,
        transport_state.spectrum_frequency.value,
        number_of_vpackets,
        montecarlo_configuration.VPACKET_LOGGING,
        iteration=iteration,
        show_progress_bars=show_progress_bars,
        total_iterations=total_iterations,
    )

    transport_state._montecarlo_virtual_luminosity.value[
        :
    ] = v_packets_energy_hist
    transport_state.last_interaction_type = last_interaction_tracker.types
    transport_state.last_interaction_in_nu = last_interaction_tracker.in_nus
    transport_state.last_line_interaction_in_id = (
        last_interaction_tracker.in_ids
    )
    transport_state.last_line_interaction_out_id = (
        last_interaction_tracker.out_ids
    )
    transport_state.last_line_interaction_shell_id = (
        last_interaction_tracker.shell_ids
    )

    if montecarlo_configuration.VPACKET_LOGGING and number_of_vpackets > 0:
        transport_state.virt_packet_nus = np.concatenate(
            virt_packet_nus
        ).ravel()
        transport_state.virt_packet_energies = np.concatenate(
            virt_packet_energies
        ).ravel()
        transport_state.virt_packet_initial_mus = np.concatenate(
            virt_packet_initial_mus
        ).ravel()
        transport_state.virt_packet_initial_rs = np.concatenate(
            virt_packet_initial_rs
        ).ravel()
        transport_state.virt_packet_last_interaction_in_nu = np.concatenate(
            virt_packet_last_interaction_in_nu
        ).ravel()
        transport_state.virt_packet_last_interaction_type = np.concatenate(
            virt_packet_last_interaction_type
        ).ravel()
        transport_state.virt_packet_last_line_interaction_in_id = (
            np.concatenate(virt_packet_last_line_interaction_in_id).ravel()
        )
        transport_state.virt_packet_last_line_interaction_out_id = (
            np.concatenate(virt_packet_last_line_interaction_out_id).ravel()
        )
        transport_state.virt_packet_last_line_interaction_shell_id = (
            np.concatenate(virt_packet_last_line_interaction_shell_id).ravel()
        )
    update_iterations_pbar(1)
    refresh_packet_pbar()
    # Condition for Checking if RPacket Tracking is enabled
    if montecarlo_configuration.RPACKET_TRACKING:
        transport_state.rpacket_tracker = rpacket_trackers


@njit(**njit_dict)
def montecarlo_main_loop(
    packet_collection,
    numba_radial_1d_geometry,
    numba_model,
    opacity_state,
    estimators,
    spectrum_frequency,
    number_of_vpackets,
    virtual_packet_logging,
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
                montecarlo_configuration.v_packet_spawn_start_frequency,
                montecarlo_configuration.v_packet_spawn_end_frequency,
                number_of_vpackets,
                montecarlo_configuration.temporary_v_packet_bins,
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
            Estimators(
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
    # Arrays for vpacket logging
    virt_packet_nus = []
    virt_packet_energies = []
    virt_packet_initial_mus = []
    virt_packet_initial_rs = []
    virt_packet_last_interaction_in_nu = []
    virt_packet_last_interaction_type = []
    virt_packet_last_line_interaction_in_id = []
    virt_packet_last_line_interaction_out_id = []
    virt_packet_last_line_interaction_shell_id = []
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

        # Get the local estimators for this thread
        vpacket_collection = vpacket_collections[i]

        # RPacket Tracker for this thread
        rpacket_tracker = rpacket_trackers[i]

        loop = single_packet_loop(
            r_packet,
            numba_radial_1d_geometry,
            numba_model,
            opacity_state,
            local_estimators,
            vpacket_collection,
            rpacket_tracker,
        )
        packet_collection.output_nus[i] = r_packet.nu
        # output_nus[i] = r_packet.nu
        last_interaction_tracker.update_last_interaction(r_packet, i)

        if r_packet.status == PacketStatus.REABSORBED:
            packet_collection.output_energies[i] = -r_packet.energy
            last_interaction_tracker.types[i] = r_packet.last_interaction_type
        elif r_packet.status == PacketStatus.EMITTED:
            packet_collection.output_energies[i] = r_packet.energy
            last_interaction_tracker.types[i] = r_packet.last_interaction_type

        vpackets_nu = vpacket_collection.nus[: vpacket_collection.idx]
        vpackets_energy = vpacket_collection.energies[: vpacket_collection.idx]
        vpackets_initial_mu = vpacket_collection.initial_mus[
            : vpacket_collection.idx
        ]
        vpackets_initial_r = vpacket_collection.initial_rs[
            : vpacket_collection.idx
        ]

        v_packets_idx = np.floor(
            (vpackets_nu - spectrum_frequency[0]) / delta_nu
        ).astype(np.int64)

        for j, idx in enumerate(v_packets_idx):
            if (vpackets_nu[j] < spectrum_frequency[0]) or (
                vpackets_nu[j] > spectrum_frequency[-1]
            ):
                continue
            v_packets_energy_hist[idx] += vpackets_energy[j]

    for sub_estimator in estimator_list:
        estimators.increment(sub_estimator)

    if virtual_packet_logging:
        for vpacket_collection in vpacket_collections:
            vpackets_nu = vpacket_collection.nus[: vpacket_collection.idx]
            vpackets_energy = vpacket_collection.energies[
                : vpacket_collection.idx
            ]
            vpackets_initial_mu = vpacket_collection.initial_mus[
                : vpacket_collection.idx
            ]
            vpackets_initial_r = vpacket_collection.initial_rs[
                : vpacket_collection.idx
            ]
            virt_packet_nus.append(np.ascontiguousarray(vpackets_nu))
            virt_packet_energies.append(np.ascontiguousarray(vpackets_energy))
            virt_packet_initial_mus.append(
                np.ascontiguousarray(vpackets_initial_mu)
            )
            virt_packet_initial_rs.append(
                np.ascontiguousarray(vpackets_initial_r)
            )
            virt_packet_last_interaction_in_nu.append(
                np.ascontiguousarray(
                    vpacket_collection.last_interaction_in_nu[
                        : vpacket_collection.idx
                    ]
                )
            )
            virt_packet_last_interaction_type.append(
                np.ascontiguousarray(
                    vpacket_collection.last_interaction_type[
                        : vpacket_collection.idx
                    ]
                )
            )
            virt_packet_last_line_interaction_in_id.append(
                np.ascontiguousarray(
                    vpacket_collection.last_interaction_in_id[
                        : vpacket_collection.idx
                    ]
                )
            )
            virt_packet_last_line_interaction_out_id.append(
                np.ascontiguousarray(
                    vpacket_collection.last_interaction_out_id[
                        : vpacket_collection.idx
                    ]
                )
            )
            virt_packet_last_line_interaction_shell_id.append(
                np.ascontiguousarray(
                    vpacket_collection.last_interaction_shell_id[
                        : vpacket_collection.idx
                    ]
                )
            )

    if montecarlo_configuration.RPACKET_TRACKING:
        for rpacket_tracker in rpacket_trackers:
            rpacket_tracker.finalize_array()

    return (
        v_packets_energy_hist,
        last_interaction_tracker,
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_initial_mus,
        virt_packet_initial_rs,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
        virt_packet_last_line_interaction_shell_id,
        rpacket_trackers,
    )
