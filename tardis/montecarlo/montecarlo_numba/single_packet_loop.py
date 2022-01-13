from numba import njit
import numpy as np

from tardis.montecarlo.montecarlo_numba.r_packet import (
    PacketStatus,
    trace_packet,
    move_packet_across_shell_boundary,
    move_r_packet,
)

from tardis.montecarlo.montecarlo_numba.utils import MonteCarloException

from tardis.montecarlo.montecarlo_numba.frame_transformations import (
    get_inverse_doppler_factor,
    get_doppler_factor
)
from tardis.montecarlo.montecarlo_numba.interaction import (
    InteractionType,
    thomson_scatter,
    line_scatter,
    continuum_event,
)
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)

from tardis.montecarlo.montecarlo_numba.vpacket import trace_vpacket_volley

from tardis import constants as const

C_SPEED_OF_LIGHT = const.c.to("cm/s").value

from tardis.io.logger.montecarlo_tracking import log_decorator
from tardis.io.logger import montecarlo_tracking as mc_tracker

# @log_decorator
@njit
def single_packet_loop(
    r_packet, numba_model, numba_plasma, estimators, vpacket_collection,
    continuum
):
    """
    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    vpacket_collection : tardis.montecarlo.montecarlo_numba.numba_interface.VPacketCollection

    Returns
    -------
    None
        This function does not return anything but changes the r_packet object
        and if virtual packets are requested - also updates the vpacket_collection
    """
    line_interaction_type = montecarlo_configuration.line_interaction_type

    if montecarlo_configuration.full_relativity:
        set_packet_props_full_relativity(r_packet, numba_model)
    else:
        set_packet_props_partial_relativity(r_packet, numba_model)
    r_packet.initialize_line_id(numba_plasma, numba_model)

    doppler_factor = get_doppler_factor(
            r_packet.r, r_packet.mu, numba_model.time_explosion
            )
    trace_vpacket_volley(
        r_packet, vpacket_collection, numba_model, numba_plasma, continuum
    )

    if mc_tracker.DEBUG_MODE:
        r_packet_track_nu = [r_packet.nu]
        r_packet_track_mu = [r_packet.mu]
        r_packet_track_r = [r_packet.r]
        r_packet_track_interaction = [InteractionType.BOUNDARY]
        r_packet_track_distance = [0.0]

    while r_packet.status == PacketStatus.IN_PROCESS:
        #print('TRACE PACKET')
        distance, interaction_type, delta_shell = trace_packet(
            r_packet, numba_model, numba_plasma,
            estimators, continuum
        )

        if interaction_type == InteractionType.BOUNDARY:
            #print("BOUNDARY")
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            move_packet_across_shell_boundary(
                r_packet, delta_shell, len(numba_model.r_inner)
            )

        elif interaction_type == InteractionType.LINE:
            #print("LINE")
            r_packet.last_interaction_type = 2
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            line_scatter(
                r_packet,
                numba_model.time_explosion,
                line_interaction_type,
                numba_plasma,
                continuum,
            )
            trace_vpacket_volley(
                r_packet, vpacket_collection, numba_model, numba_plasma, continuum
            )

        elif interaction_type == InteractionType.ESCATTERING:
            #print("ESCATTERING")
            r_packet.last_interaction_type = 1

            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            thomson_scatter(r_packet, numba_model.time_explosion)

            trace_vpacket_volley(
                r_packet, vpacket_collection, numba_model, numba_plasma, continuum
            )
            #print("Done ESCATTERING")
        elif interaction_type == InteractionType.CONTINUUM_PROCESS:
            #print("CONTINUUM_PROCESS")
            r_packet.last_interaction_type = InteractionType.CONTINUUM_PROCESS
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            continuum_event(r_packet, numba_model.time_explosion,
                    continuum, numba_plasma)

            trace_vpacket_volley(
                r_packet, vpacket_collection, numba_model, numba_plasma, continuum
            )
        else:
            #print("OTHER")
            pass


        if mc_tracker.DEBUG_MODE:
            r_packet_track_nu.append(r_packet.nu)
            r_packet_track_mu.append(r_packet.mu)
            r_packet_track_r.append(r_packet.r)
            r_packet_track_interaction.append(interaction_type)
            r_packet_track_distance.append(distance)

    if mc_tracker.DEBUG_MODE:
        return (
            r_packet_track_nu,
            r_packet_track_mu,
            r_packet_track_r,
            r_packet_track_interaction,
            r_packet_track_distance,
        )
    # check where else initialize line ID happens!


@njit
def set_packet_props_partial_relativity(r_packet, numba_model):
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        numba_model.time_explosion,
    )
    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor


@njit
def set_packet_props_full_relativity(r_packet, numba_model):
    beta = (r_packet.r / numba_model.time_explosion) / C_SPEED_OF_LIGHT

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        numba_model.time_explosion,
    )

    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor
    r_packet.mu = (r_packet.mu + beta) / (1 + beta * r_packet.mu)
