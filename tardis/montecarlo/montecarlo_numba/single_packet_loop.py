from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba.rpacket import (
    InteractionType, PacketStatus, get_doppler_factor, trace_packet,
    move_packet_across_shell_boundary, move_rpacket)
from tardis.montecarlo.montecarlo_numba.interaction import (
    general_scatter, LineInteractionType, line_scatter)

from tardis.montecarlo.montecarlo_numba.vpacket import trace_vpacket_volley




@njit
def single_packet_loop(r_packet, numba_model, numba_plasma, estimators,
                       vpacket_collection, track_rpackets=False):
    """

    Parameters
    ----------
    r_packet: tardis.montecarlo.montecarlo_numba.rpacket.RPacket
    numba_model: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators: tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    vpacket_collection: tardis.montecarlo.montecarlo_numba.numba_interface.VPacketCollection

    Returns
    -------
        : None

    This function does not return anything but changes the r_packet object
    and if virtual packets are requested - also updates the vpacket_collection

    """
    np.random.seed(r_packet.index)
    line_interaction_type = LineInteractionType.MACROATOM

    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        numba_model.time_explosion)
    r_packet.nu /= doppler_factor
    r_packet.energy /= doppler_factor
    r_packet.initialize_line_id(numba_plasma, numba_model)

    trace_vpacket_volley(r_packet, vpacket_collection, numba_model,
                         numba_plasma)

    if track_rpackets:
        rpacket_track_nu = [r_packet.nu]
        rpacket_track_mu = [r_packet.mu]
        rpacket_track_r = [r_packet.r]
        rpacket_track_interaction = [InteractionType.BOUNDARY]
        rpacket_track_distance = [0.]


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
            line_scatter(r_packet, numba_model.time_explosion,
                         line_interaction_type, numba_plasma)

            trace_vpacket_volley(
                r_packet, vpacket_collection, numba_model, numba_plasma)

        elif interaction_type == InteractionType.ESCATTERING:
            move_rpacket(r_packet, distance, numba_model.time_explosion,
                         estimators)
            general_scatter(r_packet, numba_model.time_explosion)

            trace_vpacket_volley(r_packet, vpacket_collection, numba_model,
                                 numba_plasma)

        if track_rpackets:
            rpacket_track_nu.append(r_packet.nu)
            rpacket_track_mu.append(r_packet.mu)
            rpacket_track_r.append(r_packet.r)
            rpacket_track_interaction.append(interaction_type)
            rpacket_track_distance.append(distance)


    if track_rpackets is True:
        return (rpacket_track_nu, rpacket_track_mu, rpacket_track_r,
                rpacket_track_interaction, rpacket_track_distance)

