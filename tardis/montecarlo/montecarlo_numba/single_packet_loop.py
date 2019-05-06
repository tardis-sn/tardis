from numba import njit

from tardis.montecarlo.montecarlo_numba.rpacket import (
    InteractionType, PacketStatus, get_doppler_factor, trace_packet,
    move_packet_across_shell_boundary, move_packet, scatter, initialize_line_id)

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

    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        numba_model.time_explosion)
    r_packet.nu /= doppler_factor
    r_packet.energy /= doppler_factor
    initialize_line_id(r_packet, numba_plasma, numba_model)

    while r_packet.status == PacketStatus.IN_PROCESS:
        distance, interaction_type, delta_shell = trace_packet(
            r_packet, numba_model, numba_plasma)
        if interaction_type == InteractionType.BOUNDARY:
            move_packet(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            move_packet_across_shell_boundary(r_packet, distance, delta_shell,
                                                       len(numba_model.r_inner))
        elif interaction_type == InteractionType.LINE:
            move_packet(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            scatter(r_packet, numba_model.time_explosion)
            r_packet.next_line_id += 1
        elif interaction_type == InteractionType.ESCATTERING:
            move_packet(r_packet, distance, numba_model.time_explosion,
                                 estimators)
            scatter(r_packet, numba_model.time_explosion)
