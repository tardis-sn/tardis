from numba import prange, njit, jit
import logging
import numpy as np

from tardis.montecarlo.montecarlo_numba.r_packet import (
    RPacket,
    PacketStatus,
    MonteCarloException,
)
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection,
    VPacketCollection,
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


def montecarlo_radial1d(model, plasma, runner):
    packet_collection = PacketCollection(
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
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
    ) = montecarlo_main_loop(
        packet_collection,
        numba_model,
        numba_plasma,
        estimators,
        runner.spectrum_frequency.value,
        number_of_vpackets,
        packet_seeds,
    )

    runner._montecarlo_virtual_luminosity.value[:] = v_packets_energy_hist
    runner.last_interaction_type = last_interaction_type

    if montecarlo_configuration.VPACKET_LOGGING and number_of_vpackets > 0:
        runner.virt_packet_nus = np.concatenate(
            np.array(virt_packet_nus)
        ).ravel()
        runner.virt_packet_energies = np.concatenate(
            np.array(virt_packet_energies)
        ).ravel()
        runner.virt_packet_last_interaction_in_nu = np.array(
            virt_packet_last_interaction_in_nu
        )
        runner.virt_packet_last_interaction_type = np.array(
            virt_packet_last_interaction_type
        )
        runner.virt_packet_last_line_interaction_in_id = np.array(
            virt_packet_last_line_interaction_in_id
        )
        runner.virt_packet_last_line_interaction_out_id = np.array(
            virt_packet_last_line_interaction_out_id
        )


@njit(**njit_dict, nogil=True)
def montecarlo_main_loop(
    packet_collection,
    numba_model,
    numba_plasma,
    estimators,
    spectrum_frequency,
    number_of_vpackets,
    packet_seeds,
):
    """
    This is the main loop of the MonteCarlo routine that generates packets
    and sends them through the ejecta.

    Parameters
    ----------
    packet_collection: PacketCollection
    numba_model: NumbaModel
    estimators: NumbaEstimators
    spectrum_frequency: astropy.units.Quantity
        frequency bins
    number_of_vpackets: int
        VPackets released per interaction
    packet_seeds: numpy.array
    """
    output_nus = np.empty_like(packet_collection.packets_output_nu)
    last_interaction_types = (
        np.ones_like(packet_collection.packets_output_nu) * -1
    )
    output_energies = np.empty_like(packet_collection.packets_output_nu)

    v_packets_energy_hist = np.zeros_like(spectrum_frequency)
    delta_nu = spectrum_frequency[1] - spectrum_frequency[0]

    virt_packet_nus = []
    virt_packet_energies = []
    virt_packet_last_interaction_in_nu = []
    virt_packet_last_interaction_type = []
    virt_packet_last_line_interaction_in_id = []
    virt_packet_last_line_interaction_out_id = []

    print("Running post-merge numba montecarlo (with C close lines)!")
    for i in prange(len(output_nus)):
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
            0,
        )
        vpacket_collection = VPacketCollection(
            r_packet.index,
            spectrum_frequency,
            montecarlo_configuration.v_packet_spawn_start_frequency,
            montecarlo_configuration.v_packet_spawn_end_frequency,
            number_of_vpackets,
            montecarlo_configuration.temporary_v_packet_bins,
        )
        loop = single_packet_loop(
            r_packet, numba_model, numba_plasma, estimators, vpacket_collection
        )
        # if loop and 'stop' in loop:
        #     raise MonteCarloException

        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
            last_interaction_types[i] = 0
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
            last_interaction_types[i] = 1

        vpackets_nu = vpacket_collection.nus[: vpacket_collection.idx]
        vpackets_energy = vpacket_collection.energies[: vpacket_collection.idx]

        v_packets_idx = np.floor(
            (vpackets_nu - spectrum_frequency[0]) / delta_nu
        ).astype(np.int64)
        # if we're only in a single-packet mode
        # if montecarlo_configuration.single_packet_seed == -1:
        #     break
        for j, idx in enumerate(v_packets_idx):
            if (vpackets_nu[j] < spectrum_frequency[0]) or (
                vpackets_nu[j] > spectrum_frequency[-1]
            ):
                continue
            v_packets_energy_hist[idx] += vpackets_energy[j]

        if montecarlo_configuration.VPACKET_LOGGING:
            virt_packet_nus.append(
                vpacket_collection.nus[: vpacket_collection.idx]
            )
            virt_packet_energies.append(
                vpacket_collection.energies[: vpacket_collection.idx]
            )
            virt_packet_last_interaction_in_nu.append(
                vpacket_collection.last_interaction_in_nu
            )
            virt_packet_last_interaction_type.append(
                vpacket_collection.last_interaction_type
            )
            virt_packet_last_line_interaction_in_id.append(
                vpacket_collection.last_interaction_in_id
            )
            virt_packet_last_line_interaction_out_id.append(
                vpacket_collection.last_interaction_out_id
            )

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]

    return (
        v_packets_energy_hist,
        last_interaction_types,
        virt_packet_nus,
        virt_packet_energies,
        virt_packet_last_interaction_in_nu,
        virt_packet_last_interaction_type,
        virt_packet_last_line_interaction_in_id,
        virt_packet_last_line_interaction_out_id,
    )
