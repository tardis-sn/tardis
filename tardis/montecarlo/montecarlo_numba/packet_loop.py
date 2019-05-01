from numba import njit
import numpy as np
from enum import Enum

class PacketStatus(Enum):
  IN_PROCESS = 0
  EMITTED = 1
  REABSORBED = 2

@njit
def one_packet_loop(storage, packet, virtual_packet, mt_state):
    packet.tau_event = 0.0
    packet.nu_line = 0.0
    packet.virtual_packet = virtual_packet
    packet.status = PacketStatus.IN_PROCESS
    # initializing tau_event if it's a real packet
    if virtual_packet == 0:
        packet.reset_tau_event(mt_state)
    
    # for a virtual packet tau_event is the sum of all the tau's that the packet passes 
    while packet.status == PacketStatus.IN_PROCESS:
        # Check if we are at the end of line list.
        if not packet.last_line:
            packet.nu_line = storeage.line_list_nu[packet.next_line_id]
        
        # FIXME: get_event_handler(packet, storage, &distance, mt_state) (packet, storage, distance, mt_state)
        
        if virtual_packet > 0 and packet.tau_event > storage.tau_russian:
            event_random = rk_double(mt_state)

            if event_random > storage.survival_probability:
                packet.energy = 0.0
                packet.status = PacketStatus.EMITTED
            else:
                packet.energy = packet.energy / storage.survival_probability * \
                  np.exp(-1.0 * packet.tau_event)

    if virtual_packet > 0:
        packet.energy = packet.energy * np.exp(-1.0 * packet.tau_event)

    if packet.status == PacketStatus.REABSORBED:
        return 1
    else:
        return 0



def rpacket_interactions(rpacket):
    if 
int64_t

# Void
# montecarlo_one_packet_loop (storage_model_t * storage, rpacket_t * packet,
#                             int64_t virtual_packet, rk_state *mt_state)
# {
#   rpacket_set_tau_event (packet, 0.0);
#   rpacket_set_nu_line (packet, 0.0);
#   rpacket_set_virtual_packet (packet, virtual_packet);
#   rpacket_set_status (packet, TARDIS_PACKET_STATUS_IN_PROCESS);
#   // Initializing tau_event if it's a real packet.
#   if (virtual_packet == 0)
#     {
#       rpacket_reset_tau_event (packet,mt_state);
#     }
#   // For a virtual packet tau_event is the sum of all the tau's that the packet passes.
#   while (rpacket_get_status (packet) == TARDIS_PACKET_STATUS_IN_PROCESS)
#     {
#       // Check if we are at the end of line list.
#       if (!rpacket_get_last_line (packet))
#         {
#           rpacket_set_nu_line (packet,
#                                storage->
#                                line_list_nu[rpacket_get_next_line_id
#                                (packet)]);
#         }
#       double distance;
#       get_event_handler (packet, storage, &distance, mt_state) (packet, storage,
#                                                                 distance, mt_state);
#       if (virtual_packet > 0 && rpacket_get_tau_event (packet) > storage->tau_russian)
#         {
#             double event_random = rk_double (mt_state);
#             if (event_random > storage->survival_probability)
#               {
#                 rpacket_set_energy(packet, 0.0);
#                 rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
#               }
#             else
#               {
#                 rpacket_set_energy(packet,
#                   rpacket_get_energy (packet) / storage->survival_probability *
#                   exp (-1.0 * rpacket_get_tau_event (packet)));
#                 rpacket_set_tau_event (packet, 0.0);
#               }
#           }
#     }
#   if (virtual_packet > 0)
#     {
#       rpacket_set_energy (packet,
#                           rpacket_get_energy (packet) * exp (-1.0 *
#                                                              rpacket_get_tau_event
#                                                              (packet)));
#     }
#   return rpacket_get_status (packet) ==
#     TARDIS_PACKET_STATUS_REABSORBED ? 1 : 0;
# }
