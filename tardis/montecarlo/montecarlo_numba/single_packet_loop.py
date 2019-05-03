from numba import njit
import numpy as np
from enum import Enum
from tardis.montecarlo.montecarlo_numba.rpacket import (
    ESCATTERING, BOUNDARY, LINE, IN_PROCESS, REABSORBED)

@njit
def single_packet_loop(storage_model, r_packet):
    while r_packet.status == IN_PROCESS:
        #rpacket_interactions(r_packet, storage_model)
        r_packet.initialize_line_id(storage_model)
        distance, interaction_type, delta_shell = r_packet.trace_packet(storage_model)
        if interaction_type == BOUNDARY:
            r_packet.move_packet_across_shell_boundary(distance, delta_shell, storage_model.no_of_shells)
        elif interaction_type == LINE:
            r_packet.move_packet(distance)
            r_packet.transform_energy(storage_model)
            r_packet.line_emission(storage_model)
