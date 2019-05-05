from numba import prange, njit, config
import numpy as np
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection)
from tardis.montecarlo.montecarlo_numba.storage_model import initialize_storage_model
from tardis.montecarlo.montecarlo_numba.single_packet_loop import single_packet_loop
from tardis.montecarlo.montecarlo_numba import njit_dict


def montecarlo_radial1d(model, plasma, runner):
    packet_collection = PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy
    )
    storage_model = initialize_storage_model(model, plasma, runner)
    montecarlo_main_loop(packet_collection)


@njit(**njit_dict, nogil=True)
def montecarlo_main_loop(packet_collection):
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
        r_packet = RPacket(packet_collection.r_inner[0],
                           packet_collection.packet_mus[i],
                           packet_collection.packet_nus[i],
                           packet_collection.packet_energies[i])

        single_packet_loop(storage_model, r_packet)
        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
    storage_model.output_energies[:] = output_energies[:]
    storage_model.output_nus[:] = output_nus[:]

