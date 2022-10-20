from numba import njit

from tardis.montecarlo.montecarlo_numba.r_packet import (
    PacketStatus,
)
from tardis.transport.r_packet_transport import (
    move_r_packet,
    move_packet_across_shell_boundary,
    trace_packet,
)

from tardis.transport.frame_transformations import (
    get_inverse_doppler_factor,
    get_doppler_factor,
)
from tardis.montecarlo.montecarlo_numba.interaction import (
    InteractionType,
    thomson_scatter,
    line_scatter,
    continuum_event,
)
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)

from tardis.montecarlo.montecarlo_numba.vpacket import trace_vpacket_volley

from tardis import constants as const
from tardis.montecarlo.montecarlo_numba.opacities import (
    chi_continuum_calculator,
    chi_electron_calculator,
)
from tardis.montecarlo.montecarlo_numba.estimators import (
    update_bound_free_estimators,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@njit
def single_packet_loop(
    r_packet,
    numba_radial_1d_geometry,
    numba_model,
    numba_plasma,
    estimators,
    vpacket_collection,
    rpacket_tracker,
):
    """
    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_radial_1d_geometry : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaRadial1DGeometry
    numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    vpacket_collection : tardis.montecarlo.montecarlo_numba.numba_interface.VPacketCollection
    rpacket_collection : tardis.montecarlo.montecarlo_numba.numba_interface.RPacketCollection

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

    trace_vpacket_volley(
        r_packet,
        vpacket_collection,
        numba_radial_1d_geometry,
        numba_model,
        numba_plasma,
    )

    if montecarlo_configuration.RPACKET_TRACKING:
        rpacket_tracker.track(r_packet)

    # this part of the code is temporary and will be better incorporated
    while r_packet.status == PacketStatus.IN_PROCESS:
        # Compute continuum quantities
        # trace packet (takes opacities)
        doppler_factor = get_doppler_factor(
            r_packet.r, r_packet.mu, numba_model.time_explosion
        )
        comov_nu = r_packet.nu * doppler_factor
        chi_e = chi_electron_calculator(
            numba_plasma, comov_nu, r_packet.current_shell_id
        )
        if montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED:
            (
                chi_bf_tot,
                chi_bf_contributions,
                current_continua,
                x_sect_bfs,
                chi_ff,
            ) = chi_continuum_calculator(
                numba_plasma, comov_nu, r_packet.current_shell_id
            )
            chi_continuum = chi_e + chi_bf_tot + chi_ff
            escat_prob = chi_e / chi_continuum  # probability of e-scatter
            distance, interaction_type, delta_shell = trace_packet(
                r_packet,
                numba_radial_1d_geometry,
                numba_model,
                numba_plasma,
                estimators,
                chi_continuum,
                escat_prob,
            )
            update_bound_free_estimators(
                comov_nu,
                r_packet.energy * doppler_factor,
                r_packet.current_shell_id,
                distance,
                estimators,
                numba_plasma.t_electrons[r_packet.current_shell_id],
                x_sect_bfs,
                current_continua,
                numba_plasma.bf_threshold_list_nu,
            )
        else:
            escat_prob = 1.0
            chi_continuum = chi_e
            distance, interaction_type, delta_shell = trace_packet(
                r_packet,
                numba_radial_1d_geometry,
                numba_model,
                numba_plasma,
                estimators,
                chi_continuum,
                escat_prob,
            )

        # If continuum processes: update continuum estimators

        if interaction_type == InteractionType.BOUNDARY:
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            move_packet_across_shell_boundary(
                r_packet, delta_shell, len(numba_radial_1d_geometry.r_inner)
            )

        elif interaction_type == InteractionType.LINE:
            r_packet.last_interaction_type = 2
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            line_scatter(
                r_packet,
                numba_model.time_explosion,
                line_interaction_type,
                numba_plasma,
            )
            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                numba_plasma,
            )

        elif interaction_type == InteractionType.ESCATTERING:
            r_packet.last_interaction_type = 1

            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            thomson_scatter(r_packet, numba_model.time_explosion)

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                numba_plasma,
            )
        elif (
            montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
            and interaction_type == InteractionType.CONTINUUM_PROCESS
        ):
            r_packet.last_interaction_type = InteractionType.CONTINUUM_PROCESS
            move_r_packet(
                r_packet, distance, numba_model.time_explosion, estimators
            )
            continuum_event(
                r_packet,
                numba_model.time_explosion,
                numba_plasma,
                chi_bf_tot,
                chi_ff,
                chi_bf_contributions,
                current_continua,
            )

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                numba_plasma,
            )
        else:
            pass
        if montecarlo_configuration.RPACKET_TRACKING:
            rpacket_tracker.track(r_packet)


@njit
def set_packet_props_partial_relativity(r_packet, numba_model):
    """Sets properties of the packets given partial relativity

    Parameters
    ----------
        r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
        numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    """
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        numba_model.time_explosion,
    )
    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor


@njit
def set_packet_props_full_relativity(r_packet, numba_model):
    """Sets properties of the packets given full relativity

    Parameters
    ----------
        r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
        numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    """
    beta = (r_packet.r / numba_model.time_explosion) / C_SPEED_OF_LIGHT

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        numba_model.time_explosion,
    )

    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor
    r_packet.mu = (r_packet.mu + beta) / (1 + beta * r_packet.mu)
