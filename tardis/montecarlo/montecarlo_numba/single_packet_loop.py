from numba import njit
import numpy as np
from enum import Enum
from tardis.montecarlo.montecarlo_numba.rpacket import (IN_PROCESS, REABSORBED)
from tardis.montecarlo.montecarlo_numba.rpacket import InteractionType

@njit
def single_packet_loop(storage_model, r_packet):
    r_packet.initialize_line_id(storage_model)

    while r_packet.status == IN_PROCESS:
        distance, interaction_type, delta_shell = r_packet.trace_packet(storage_model)
        if interaction_type == InteractionType.BOUNDARY:
            r_packet.move_packet_across_shell_boundary(distance, delta_shell, storage_model.no_of_shells)
        elif interaction_type == InteractionType.LINE:
            r_packet.move_packet(distance)
            r_packet.transform_energy(storage_model)
            r_packet.line_emission(storage_model)
