from numba import prange, njit, jit
import logging
import numpy as np

from tardis.montecarlo.montecarlo_numba.r_packet import (
    RPacket,
    PacketStatus,
)
from tardis.montecarlo.montecarlo_numba.utils import MonteCarloException

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
from numba.typed import List


def montecarlo_radial1d(model, plasma, runner):
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
    runner.last_interaction_in_nu = last_interaction_in_nu
    runner.last_line_interaction_in_id = last_line_interaction_in_id
    runner.last_line_interaction_out_id = last_line_interaction_out_id

    if montecarlo_configuration.VPACKET_LOGGING and number_of_vpackets > 0:
        runner.virt_packet_nus = np.concatenate(
            np.array(virt_packet_nus, dtype=object)
        ).ravel()
        runner.virt_packet_energies = np.concatenate(
            np.array(virt_packet_energies, dtype=object)
        ).ravel()
        runner.virt_packet_initial_mus = np.concatenate(
            np.array(virt_packet_initial_mus, dtype=object)
        ).ravel()
        runner.virt_packet_initial_rs = np.concatenate(
            np.array(virt_packet_initial_rs, dtype=object)
        ).ravel()
        runner.virt_packet_last_interaction_in_nu = np.concatenate(
            np.array(virt_packet_last_interaction_in_nu, dtype=object)
        ).ravel()
        runner.virt_packet_last_interaction_type = np.concatenate(
            np.array(virt_packet_last_interaction_type, dtype=object)
        ).ravel()
        runner.virt_packet_last_line_interaction_in_id = np.concatenate(
            np.array(virt_packet_last_line_interaction_in_id, dtype=object)
        ).ravel()
        runner.virt_packet_last_line_interaction_out_id = np.concatenate(
            np.array(virt_packet_last_line_interaction_out_id, dtype=object)
        ).ravel()


@njit(**njit_dict)
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
    packet_collection : PacketCollection
    numba_model : NumbaModel
    estimators : NumbaEstimators
    spectrum_frequency : astropy.units.Quantity
        frequency bins
    number_of_vpackets : int
        VPackets released per interaction
    packet_seeds : numpy.array
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
        )
        vpacket_collection = vpacket_collections[i]

        loop = single_packet_loop(
            r_packet, numba_model, numba_plasma, estimators, vpacket_collection
        )
        # if loop and 'stop' in loop:
        #     raise MonteCarloException

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

    if montecarlo_configuration.VPACKET_LOGGING:
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
    )
