from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba.rpacket import (
    InteractionType, PacketStatus)

@njit
def single_packet_loop(r_packet, numba_model, numba_plasma, estimators):
    """

    Parameters
    ----------
    r_packet: tardis.montecarlo.montecarlo_numba.rpacket.RPacket
    numba_model: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators: tardis.montecarlo.montecarlo_numba.numba_interface.Estimators

    Returns
    -------

    """
    r_packet.initialize_line_id(numba_plasma, numba_model)

    while r_packet.status == PacketStatus.IN_PROCESS:
        distance, interaction_type, delta_shell = r_packet.trace_packet(
            numba_model, numba_plasma)
        if interaction_type == InteractionType.BOUNDARY:
            r_packet.move_packet(distance, numba_model.time_explosion,
                                 estimators)
            r_packet.move_packet_across_shell_boundary(distance, delta_shell,
                                                       len(numba_model.r_inner))
        elif interaction_type == InteractionType.LINE:
            r_packet.move_packet(distance, numba_model.time_explosion,
                                 estimators)
            r_packet.scatter(numba_model.time_explosion)
            r_packet.next_line_id += 1
        elif interaction_type == InteractionType.ESCATTERING:
            r_packet.move_packet(distance, numba_model.time_explosion,
                                 estimators)
            r_packet.scatter(numba_model.time_explosion)
