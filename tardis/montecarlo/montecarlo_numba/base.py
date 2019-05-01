
from numba import prange, njit
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.storage_model import StorageModel
from tardis.montecarlo.montecarlo_numba.packet_loop import one_packet_loop



def montecarlo_radial1d(model, plasma, runner):
    storage_model = initialize_storage_model(model, plasma, runner)
    montecarlo_main_loop(storage_model)
    
    pass
        #montecarlo.montecarlo_radial1d(
        #    model, plasma, self,
        #    virtual_packet_flag=no_of_virtual_packets,
        #    nthreads=nthreads,
        #    last_run=last_run)
        # Workaround so that j_blue_estimator is in the right ordering
        # They are written as an array of dimension (no_of_shells, no_of_lines)
        # but python expects (no_of_lines, no_of_shells)

@njit
def montecarlo_main_loop(storage_model):
    for i in prange(storage.no_of_packets):
        r_packet = RPacket(r=storage.r_inner[0],
                           mu=storage.packet_nus[i],
                           nu=storage.packet_mus[i],
                           energy=storage.packet_energies[i])
        r_packet_final_state = one_packet_loop(r_packet)
        storage_model.output_nus[i] = rpacket.nu
        if r_packet_final_state == PacketStatus.REABSORBED:
            storage_model.output_energies[i] = -r_packet.energy
        elif r_packet_final_state == PacketStatus.EMITTED:
            storage_model.output_energies[i] = r_packet.energy
    return None
