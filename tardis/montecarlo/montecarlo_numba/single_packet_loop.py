from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba.rpacket import (
    InteractionType, PacketStatus)

@njit
def single_packet_loop(storage_model, r_packet):
    r_packet.initialize_line_id(storage_model)

    while r_packet.status == PacketStatus.IN_PROCESS:
        distance, interaction_type, delta_shell = r_packet.trace_packet(storage_model)
        if interaction_type == InteractionType.BOUNDARY:
            r_packet.move_packet_across_shell_boundary(distance, delta_shell, storage_model.no_of_shells)
        elif interaction_type == InteractionType.LINE:
            r_packet.move_packet(distance)
            r_packet.scatter(storage_model)
            r_packet.next_line_id += 1
        elif interaction_type == InteractionType.ESCATTERING:
            r_packet.move_packet(distance)
            r_packet.scatter(storage_model)
