from numba import prange, njit, config, int64
import numpy as np
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection, VPacketCollection, NumbaModel, numba_plasma_initialize, Estimators, MonteCarloConfiguration)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop)
from tardis.montecarlo.montecarlo_numba import njit_dict


def montecarlo_radial1d(model, plasma, runner, no_of_virtual_packets):
    packet_collection = PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy
    )
    montecarlo_configuration = MonteCarloConfiguration(no_of_virtual_packets)
    numba_model = NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)
    numba_plasma = numba_plasma_initialize(plasma)
    estimators = Estimators(runner.j_estimator, runner.nu_bar_estimator)

    v_packets_energy_hist = montecarlo_main_loop(packet_collection, numba_model, numba_plasma, estimators, runner.spectrum_frequency.value, montecarlo_configuration)
    
    runner._montecarlo_virtual_luminosity.value[:] = v_packets_energy_hist


@njit(**njit_dict, nogil=True)
def montecarlo_main_loop(packet_collection, numba_model, numba_plasma,
                         estimators, spectrum_frequency, montecarlo_configuration):
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
                           packet_collection.packets_input_energy[i])
        
        vpacket_collection = VPacketCollection(spectrum_frequency, montecarlo_configuration.number_of_vpackets, 20000)
        single_packet_loop(r_packet, numba_model, numba_plasma, estimators, vpacket_collection)

        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
        
        vpackets_nu = vpacket_collection.nus[:vpacket_collection.idx]
        vpackets_energy = vpacket_collection.energies[:vpacket_collection.idx]

        v_packets_idx = np.floor((vpackets_nu - spectrum_frequency[0]) / delta_nu).astype(np.int64)
        for i, idx in enumerate(v_packets_idx):
            v_packets_energy_hist[idx] += vpackets_energy[i]

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]
    
    return v_packets_energy_hist

