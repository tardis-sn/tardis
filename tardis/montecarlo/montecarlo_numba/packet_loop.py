from numba import njit
import numpy as np
from enum import Enum
from tardis.montecarlo.montecarlo_numba.rpacket import (
    ESCATTERING, BOUNDARY, LINE, IN_PROCESS, REABSORBED)
@njit
def one_packet_loop(storage_model, r_packet):
    #r_packet.nu_line = 0.0
    #packet.virtual_packet = virtual_packet
    r_packet.status = IN_PROCESS
    # initializing tau_event if it's a real packet
    
    # for a virtual packet tau_event is the sum of all the tau's that the packet passes 
    while r_packet.status == IN_PROCESS:
        rpacket_interactions(r_packet, storage_model)
        # Check if we are at the end of line list.
        #if not packet.last_line:
        #    packet.nu_line = storeage.line_list_nu[packet.next_line_id]
        
        # FIXME: get_event_handler(packet, storage, &distance, mt_state) (packet, storage, distance, mt_state)
        
    #    if virtual_packet > 0 and packet.tau_event > storage.tau_russian:
    #        event_random = rk_double(mt_state)
    #        if event_random > storage.survival_probability:
    #            packet.energy = 0.0
    #            packet.status = PacketStatus.EMITTED
    #        else:
    #            packet.energy = packet.energy / storage.survival_probability * \
    #              np.exp(-1.0 * packet.tau_event)
    #if virtual_packet > 0:
    #    packet.energy = packet.energy * np.exp(-1.0 * packet.tau_event)
@njit
def rpacket_interactions(r_packet, storage_model):
    r_packet.compute_distances(storage_model)
    if r_packet.next_interaction == BOUNDARY:
        r_packet.move_packet_across_shell_boundary(storage_model)
    else:
        r_packet.line_scatter(storage_model)
        
    #    pass
    #else:
    #    pass