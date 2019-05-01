from numba import prange, njit
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, REABSORBED, EMITTED
from tardis.montecarlo.montecarlo_numba.storage_model import StorageModel, initialize_storage_model
from tardis.montecarlo.montecarlo_numba.packet_loop import one_packet_loop



def montecarlo_radial1d(model, plasma, runner):
    storage_model = initialize_storage_model(model, plasma, runner)
    montecarlo_main_loop(storage_model)

@njit
def montecarlo_main_loop(storage_model):

    for i in prange(storage_model.no_of_packets):
        r_packet = RPacket(storage_model.r_inner[0],
                           storage_model.packet_mus[i],
                           storage_model.packet_nus[i],
                           storage_model.packet_energies[i])
        r_packet.compute_distances(storage_model)
        one_packet_loop(storage_model, r_packet)
        storage_model.output_nus[i] = r_packet.nu
        #print('!!!! THIS IS THE FINAL STATE !!!!!!!', r_packet_final_state)
        #return
        if r_packet.status == REABSORBED:
            storage_model.output_energies[i] = -r_packet.energy
        elif r_packet.status == EMITTED:
            storage_model.output_energies[i] = r_packet.energy
