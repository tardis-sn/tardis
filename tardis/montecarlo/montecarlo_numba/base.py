from numba import prange, njit, jit
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from tardis.montecarlo.montecarlo_numba.r_packet import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection, VPacketCollection, NumbaModel, numba_plasma_initialize,
    Estimators, MonteCarloConfiguration, configuration_initialize)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop)
from tardis.montecarlo.montecarlo_numba import njit_dict


def montecarlo_radial1d(model, plasma, runner, montecarlo_configuration):
    packet_collection = PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy
    )

    numba_model = NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)
    numba_plasma = numba_plasma_initialize(plasma)
    estimators = Estimators(runner.j_estimator, runner.nu_bar_estimator,
                            runner.j_b_lu_estimator, runner.edot_lu_estimator)

    v_packets_energy_hist = montecarlo_main_loop(
        packet_collection, numba_model, numba_plasma, estimators,
        runner.spectrum_frequency.value, montecarlo_configuration)
    
    runner._montecarlo_virtual_luminosity.value[:] = v_packets_energy_hist


@jit(**njit_dict, nogil=True)
def montecarlo_main_loop(packet_collection, numba_model, numba_plasma,
                         estimators, spectrum_frequency,
                         montecarlo_configuration):
    """
    This is the main loop of the MonteCarlo routine that generates packets 
    and sends them through the ejecta. 

    Parameters
    ----------
    storage_model : [type]
        [description]
    """

    output_nus = np.empty_like(packet_collection.packets_output_nu)
    output_energies = np.empty_like(packet_collection.packets_output_nu)

    v_packets_energy_hist = np.zeros_like(spectrum_frequency)
    delta_nu = spectrum_frequency[1] - spectrum_frequency[0]

    for i in prange(len(output_nus)):
        r_packet = RPacket(numba_model.r_inner[0],
                           packet_collection.packets_input_mu[i],
                           packet_collection.packets_input_nu[i],
                           packet_collection.packets_input_energy[i],
                           i)
        np.random.seed(r_packet.index)
        vpacket_collection = VPacketCollection(
            spectrum_frequency, montecarlo_configuration.number_of_vpackets,
            montecarlo_configuration.temporary_v_packet_bins)
        single_packet_loop(r_packet, numba_model, numba_plasma, estimators,
                           vpacket_collection,
                           montecarlo_configuration)

        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
        
        vpackets_nu = vpacket_collection.nus[:vpacket_collection.idx]
        vpackets_energy = vpacket_collection.energies[:vpacket_collection.idx]

        v_packets_idx = np.floor((vpackets_nu - spectrum_frequency[0]) /
                                 delta_nu).astype(np.int64)
        for j, idx in enumerate(v_packets_idx):
            if ((vpackets_nu[j] < spectrum_frequency[0]) or
                    (vpackets_nu[j] > spectrum_frequency[-1])):
                continue
            v_packets_energy_hist[idx] += vpackets_energy[j]

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]
    
    return v_packets_energy_hist

def plot_single_packet(runner, model, plasma, no_of_virtual_packets):
    """
    Used to track the progress of single packets.
    :param r_track_r: numpy array. the r of the tracked packet
    :param r_track_distance: the traveled distance fo the tracked packet.
    :return:
    """
    def follow_single_rpacket(i):
        r_packet = RPacket(model.r_inner[0].value,
                           packet_collection.packets_input_mu[i],
                           packet_collection.packets_input_nu[i],
                           packet_collection.packets_input_energy[i],
                           i)
        np.random.seed(i)
        r_track_nu, r_track_mu, r_track_r, r_track_interaction, r_track_distance = single_packet_loop(
            r_packet, numba_model, numba_plasma, estimators, vpacket_collection,
            montecarlo_configuration, track_rpackets=True)
        r_track_r = np.array(r_track_r)
        r_track_distance = np.array(r_track_distance)
        angle = np.cumsum(np.arccos((r_track_r[:-1] ** 2 + r_track_r[
                                                           1:] ** 2 - r_track_distance[
                                                                      1:] ** 2) / (
                                                2 * r_track_r[:-1] * r_track_r[
                                                                     1:])))
        x, y = np.sin(angle) * r_track_r[1:], np.cos(angle) * r_track_r[1:]
        x = np.hstack(([0], x))
        y = np.hstack((r_track_r[0], y))
        return x, y


    packet_collection = PacketCollection(runner.input_nu, runner.input_mu,
                                         runner.input_energy,
                                         runner._output_nu,
                                         runner._output_energy)
    numba_plasma = numba_plasma_initialize(plasma)
    numba_model = NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)
    estimators = Estimators(runner.j_estimator, runner.nu_bar_estimator,
                            runner.j_b_lu_estimator, runner.edot_lu_estimator)
    spectrum_frequency = runner.spectrum_frequency.value
    vpacket_collection = VPacketCollection(
        spectrum_frequency=spectrum_frequency,
        number_of_vpackets=2, temporary_v_packet_bins=20000)
    montecarlo_configuration = configuration_initialize(runner,
                                                        no_of_virtual_packets,
                                                        full_relativity=True)

    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    ax.set_ylim(0, 1.1 * model.radius[-1].cgs.value)
    ax.set_xlim(0, 1.1 * model.radius[-1].cgs.value)
    for r in model.radius:
        print(r.value)
        ax.add_patch(mpl.patches.Circle((0, 0), r.value, edgecolor='black',
                                        facecolor='none'))

    for packet_idx in np.random.randint(10000, size=20):
        try:
            x, y = follow_single_rpacket(packet_idx)
            ax.plot([x[0]], [y[0]], marker='o', color='red')
            ax.plot(x, y, marker='x', color='blue')
        except:
            plt.savefig('tardis_viz.pdf')
            plt.show()
            break

