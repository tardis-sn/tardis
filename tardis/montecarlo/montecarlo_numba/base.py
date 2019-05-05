from numba import prange, njit, config
import numpy as np
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection, NumbaModel, NumbaPlasma)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop)
from tardis.montecarlo.montecarlo_numba import njit_dict


def montecarlo_radial1d(model, plasma, runner):
    packet_collection = PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy
    )
    numba_model = NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)
    numba_plasma = NumbaPlasma(plasma.electron_densities.values,
                               plasma.atomic_data.lines.nu.values,
                               np.ascontiguousarray(
                                   plasma.tau_sobolevs.values.copy(),
                                   dtype=np.float64)
                               )

    #storage_model = initialize_storage_model(model, plasma, runner)
    montecarlo_main_loop(packet_collection, numba_model, numba_plasma)


@njit(**njit_dict, nogil=True)
def montecarlo_main_loop(packet_collection, numba_model, numba_plasma):
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

    for i in prange(len(output_nus)):
        r_packet = RPacket(numba_model.r_inner[0],
                           packet_collection.packets_input_mu[i],
                           packet_collection.packets_input_nu[i],
                           packet_collection.packets_input_energy[i])

        single_packet_loop(r_packet, numba_model, numba_plasma)
        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]

