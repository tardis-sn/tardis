import logging
import numpy as np
import pandas as pd

from tardis.montecarlo.montecarlo_numba.r_packet import (
    RPacket,
    PacketStatus,
)
from tardis.montecarlo.montecarlo_numba.utils import MonteCarloException

from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection,
    VPacketCollection,
    RPacketTracker,
    NumbaModel,
    numba_plasma_initialize,
    Estimators,
    configuration_initialize,
)

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop,
)
from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis.util.base import update_iterations_pbar, update_packet_pbar

from numba import prange, njit, jit, objmode
from numba.typed import List
from itertools import chain


def montecarlo_radial1d(
    model,
    plasma,
    iteration,
    no_of_packets,
    total_iterations,
    show_progress_bars,
    runner,
):
    packet_collection = PacketCollection(
        runner.input_r,
        runner.input_nu,
        runner.input_mu,
        runner.input_energy,
        runner._output_nu,
        runner._output_energy,
    )

    numba_model = NumbaModel(
        runner.r_inner_cgs,
        runner.r_outer_cgs,
        model.time_explosion.to("s").value,
    )
    numba_plasma = numba_plasma_initialize(plasma, runner.line_interaction_type)
    estimators = Estimators(
        runner.j_estimator,
        runner.nu_bar_estimator,
        runner.j_blue_estimator,
        runner.Edotlu_estimator,
    )

    packet_seeds = montecarlo_configuration.packet_seeds

    number_of_vpackets = montecarlo_configuration.number_of_vpackets

    (
        v_packets_energy_hist,
        last_interaction_type,
        last_interaction_in_nu,
        last_line_interaction_in_id,
        last_line_interaction_out_id,
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_initial_mus,
        virt_packet_initial_rs,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
        rpacket_collections,
    ) = montecarlo_main_loop(
        packet_collection,
        numba_model,
        numba_plasma,
        estimators,
        runner.spectrum_frequency.value,
        number_of_vpackets,
        packet_seeds,
        montecarlo_configuration.VPACKET_LOGGING,
        iteration=iteration,
        show_progress_bars=show_progress_bars,
        no_of_packets=no_of_packets,
        total_iterations=total_iterations,
    )

    runner._montecarlo_virtual_luminosity.value[:] = v_packets_energy_hist
    runner.last_interaction_type = last_interaction_type
    runner.last_interaction_in_nu = last_interaction_in_nu
    runner.last_line_interaction_in_id = last_line_interaction_in_id
    runner.last_line_interaction_out_id = last_line_interaction_out_id

    if montecarlo_configuration.VPACKET_LOGGING and number_of_vpackets > 0:
        runner.virt_packet_nus = np.concatenate(
            np.array(virt_packet_nus)
        ).ravel()
        runner.virt_packet_energies = np.concatenate(
            np.array(virt_packet_energies)
        ).ravel()
        runner.virt_packet_initial_mus = np.concatenate(
            np.array(virt_packet_initial_mus)
        ).ravel()
        runner.virt_packet_initial_rs = np.concatenate(
            np.array(virt_packet_initial_rs)
        ).ravel()
        runner.virt_packet_last_interaction_in_nu = np.concatenate(
            np.array(virt_packet_last_interaction_in_nu)
        ).ravel()
        runner.virt_packet_last_interaction_type = np.concatenate(
            np.array(virt_packet_last_interaction_type)
        ).ravel()
        runner.virt_packet_last_line_interaction_in_id = np.concatenate(
            np.array(virt_packet_last_line_interaction_in_id)
        ).ravel()
        runner.virt_packet_last_line_interaction_out_id = np.concatenate(
            np.array(virt_packet_last_line_interaction_out_id)
        ).ravel()
    update_iterations_pbar(1)

    # Condition for Checking if R Packet Tracking is enabled
    if montecarlo_configuration.RPACKET_TRACKING:
        # Creates a dataframe for a particular rpacket
        runner.rpacket_collections = rpacket_collections
        tracked_df = create_tracked_rpacket_df(
            rpacket_collections,
            iteration,
        )

        # Appending the dataframe based on iterations
        runner.rpacket_tracker = track_rpacket_dataframe(
            runner.rpacket_tracker, tracked_df
        )


@njit(**njit_dict)
def montecarlo_main_loop(
    packet_collection,
    numba_model,
    numba_plasma,
    estimators,
    spectrum_frequency,
    number_of_vpackets,
    packet_seeds,
    virtual_packet_logging,
    iteration,
    show_progress_bars,
    no_of_packets,
    total_iterations,
):
    """
    This is the main loop of the MonteCarlo routine that generates packets
    and sends them through the ejecta.

    Parameters
    ----------
    packet_collection : PacketCollection
    numba_model : NumbaModel
    estimators : NumbaEstimators
    spectrum_frequency : astropy.units.Quantity
        frequency bins
    number_of_vpackets : int
        VPackets released per interaction
    packet_seeds : numpy.array
    virtual_packet_logging : bool
        Option to enable virtual packet logging.
    """
    output_nus = np.empty_like(packet_collection.packets_output_nu)
    last_interaction_types = (
        np.ones_like(packet_collection.packets_output_nu, dtype=np.int64) * -1
    )
    output_energies = np.empty_like(packet_collection.packets_output_nu)

    last_interaction_in_nus = np.empty_like(packet_collection.packets_output_nu)
    last_line_interaction_in_ids = (
        np.ones_like(packet_collection.packets_output_nu, dtype=np.int64) * -1
    )
    last_line_interaction_out_ids = (
        np.ones_like(packet_collection.packets_output_nu, dtype=np.int64) * -1
    )

    v_packets_energy_hist = np.zeros_like(spectrum_frequency)
    delta_nu = spectrum_frequency[1] - spectrum_frequency[0]

    # Pre-allocate a list of vpacket collections for later storage
    vpacket_collections = List()
    for i in range(len(output_nus)):
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

    # Arrays for vpacket logging
    virt_packet_nus = []
    virt_packet_energies = []
    virt_packet_initial_mus = []
    virt_packet_initial_rs = []
    virt_packet_last_interaction_in_nu = []
    virt_packet_last_interaction_type = []
    virt_packet_last_line_interaction_in_id = []
    virt_packet_last_line_interaction_out_id = []

    # Configuring the Tracking for R_Packets
    rpacket_collections = List()
    for i in range(len(output_nus)):
        rpacket_collections.append(RPacketTracker())

    for i in prange(len(output_nus)):
        if show_progress_bars:
            with objmode:
                update_amount = 1
                update_packet_pbar(
                    update_amount,
                    current_iteration=iteration,
                    no_of_packets=no_of_packets,
                    total_iterations=total_iterations,
                )

        if montecarlo_configuration.single_packet_seed != -1:
            seed = packet_seeds[montecarlo_configuration.single_packet_seed]
            np.random.seed(seed)
        else:
            seed = packet_seeds[i]
            np.random.seed(seed)
        r_packet = RPacket(
            numba_model.r_inner[0],
            packet_collection.packets_input_mu[i],
            packet_collection.packets_input_nu[i],
            packet_collection.packets_input_energy[i],
            seed,
            i,
        )

        vpacket_collection = vpacket_collections[i]
        rpacket_collection = rpacket_collections[i]

        single_packet_loop(
            r_packet,
            numba_model,
            numba_plasma,
            estimators,
            vpacket_collection,
            rpacket_collection,
        )

        output_nus[i] = r_packet.nu
        last_interaction_in_nus[i] = r_packet.last_interaction_in_nu
        last_line_interaction_in_ids[i] = r_packet.last_line_interaction_in_id
        last_line_interaction_out_ids[i] = r_packet.last_line_interaction_out_id

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
            last_interaction_types[i] = r_packet.last_interaction_type
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
            last_interaction_types[i] = r_packet.last_interaction_type

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
        # if we're only in a single-packet mode
        # if montecarlo_configuration.single_packet_seed == -1:
        #    break
        for j, idx in enumerate(v_packets_idx):
            if (vpackets_nu[j] < spectrum_frequency[0]) or (
                vpackets_nu[j] > spectrum_frequency[-1]
            ):
                continue
            v_packets_energy_hist[idx] += vpackets_energy[j]

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

    if montecarlo_configuration.RPACKET_TRACKING:
        for i in range(len(rpacket_collections)):
            rpacket_collections[i].finalize_array()

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]

    return (
        v_packets_energy_hist,
        last_interaction_types,
        last_interaction_in_nus,
        last_line_interaction_in_ids,
        last_line_interaction_out_ids,
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_initial_mus,
        virt_packet_initial_rs,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
        rpacket_collections,
    )


def create_tracked_rpacket_df(rpacket_collections, iteration):
    """
    Creates Dataframe for particular rpacket along with Iteration info

    Parameters
    ----------
    rpacket_collections : list (contains `numba.jitclass` instances)
        A list of rpackets. Stores all the data for each interaction a particular rpacket undergoes as the simulation progresses.
    iteration : int
        Current Simulation 'Iteration' value

    Returns
    -------
    rpacket_tracked_df : pd.DataFrame
        A dataframe that is generated with all the properties of a single rpacket.
    """

    seed = rpacket_collections[0].seed
    index = rpacket_collections[0].index
    status = rpacket_collections[0].status
    r = rpacket_collections[0].r
    nu = rpacket_collections[0].nu
    mu = rpacket_collections[0].mu
    energy = rpacket_collections[0].energy
    shell_id = rpacket_collections[0].shell_id

    for i in range(1, len(rpacket_collections)):
        seed = chain(seed, rpacket_collections[i].seed)
        index = chain(index, rpacket_collections[i].index)
        status = chain(status, rpacket_collections[i].status)
        r = chain(r, rpacket_collections[i].r)
        nu = chain(nu, rpacket_collections[i].nu)
        mu = chain(mu, rpacket_collections[i].mu)
        energy = chain(energy, rpacket_collections[i].energy)
        shell_id = chain(shell_id, rpacket_collections[i].shell_id)

    vals = [
        index,
        seed,
        status,
        r,
        nu,
        mu,
        energy,
        shell_id,
    ]
    columns_name = [
        "Packet Index",
        "Packet Seed",
        "Packet Status",
        "r",
        "nu",
        "mu",
        "energy",
        "shell_id",
    ]

    rpacket_tracked_df = pd.DataFrame(
        zip(*vals), columns=columns_name, dtype=object
    )

    rpacket_tracked_df["Iteration"] = iteration
    return rpacket_tracked_df


def track_rpacket_dataframe(tracked_rpacket_df, tracked_df):
    """
    Appending function for joining the pre-exisiting dataframe stored in `runner.rpacket_tracker` with the newly created dataframe with iteration information.
    Helps create the final dataframe which has all the data stored with iteration for all the rpackets.

    Parameter
    ---------
    tracked_rpacket_df : pd.DataFrame
        The final DataFrame that is stores all the data for all the interactions for all the rpackets, iteration wise.
    tracked_df : pd.DataFrame
        Particular iteraction based DataFrame that stores the properties of a particular `rpacket` under consideration.

    Returns
    -------
    tracked_rpacket_df : pd.DataFrame
        Returns the Final DataFrame after appending the data of `tracked_df` into the pre-exisiting values.
    """
    tracked_rpacket_df = tracked_rpacket_df.append(
        tracked_df, ignore_index=True
    )
    tracked_rpacket_df = tracked_rpacket_df.convert_dtypes()
    columns_reorder = [
        "Iteration",
        "Packet Index",
        "Packet Seed",
        "Packet Status",
        "r",
        "nu",
        "mu",
        "energy",
        "shell_id",
    ]
    tracked_rpacket_df = tracked_rpacket_df[columns_reorder]
    return tracked_rpacket_df
