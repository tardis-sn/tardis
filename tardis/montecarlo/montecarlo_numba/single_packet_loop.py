from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba.rpacket import (
    InteractionType, PacketStatus, get_doppler_factor, trace_packet,
    move_packet_across_shell_boundary, move_rpacket, scatter)

from tardis.montecarlo.montecarlo_numba.vpacket import trace_vpacket_volley
from tardis.montecarlo.montecarlo_numba.numba_interface import VPacketCollection
@njit
def single_packet_loop(r_packet, numba_model, numba_plasma, estimators, vpacket_collection):
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

    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        numba_model.time_explosion)
    r_packet.nu /= doppler_factor
    r_packet.energy /= doppler_factor
    r_packet.initialize_line_id(numba_plasma, numba_model)

    trace_vpacket_volley(r_packet, vpacket_collection, numba_model, numba_plasma)

    while r_packet.status == PacketStatus.IN_PROCESS:
        distance, interaction_type, delta_shell = trace_packet(
            r_packet, numba_model, numba_plasma)
        if interaction_type == InteractionType.BOUNDARY:
            move_rpacket(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            move_packet_across_shell_boundary(r_packet, delta_shell,
                                                       len(numba_model.r_inner))
        elif interaction_type == InteractionType.LINE:
            move_rpacket(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            scatter(r_packet, numba_model.time_explosion)
            r_packet.next_line_id += 1
            
            trace_vpacket_volley(r_packet, vpacket_collection, numba_model, numba_plasma)

        elif interaction_type == InteractionType.ESCATTERING:
            move_rpacket(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            scatter(r_packet, numba_model.time_explosion)

            trace_vpacket_volley(r_packet, vpacket_collection, numba_model, numba_plasma)

