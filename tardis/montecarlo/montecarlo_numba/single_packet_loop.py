from numba import njit

from tardis import constants as const
from tardis.montecarlo.estimators.radfield_estimator_calcs import (
    update_bound_free_estimators,
)
from tardis.montecarlo.montecarlo_numba.interaction import (
    continuum_event,
    line_scatter,
    thomson_scatter,
)
from tardis.montecarlo.montecarlo_numba.opacities import (
    chi_continuum_calculator,
    chi_electron_calculator,
)
from tardis.montecarlo.montecarlo_numba.r_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.montecarlo.montecarlo_numba.vpacket import trace_vpacket_volley
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.r_packet_transport import (
    move_packet_across_shell_boundary,
    move_r_packet,
    trace_packet,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@njit
def single_packet_loop(
    r_packet,
    numba_radial_1d_geometry,
    numba_model,
    opacity_state,
    estimators,
    vpacket_collection,
    rpacket_tracker,
    montecarlo_configuration,
):
    """
    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_radial_1d_geometry : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaRadial1DGeometry
    numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    opacity_state : tardis.montecarlo.montecarlo_numba.numba_interface.OpacityState
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    vpacket_collection : tardis.montecarlo.montecarlo_numba.numba_interface.VPacketCollection
    rpacket_collection : tardis.montecarlo.montecarlo_numba.numba_interface.RPacketCollection

    Returns
    -------
    None
        This function does not return anything but changes the r_packet object
        and if virtual packets are requested - also updates the vpacket_collection
    """
    line_interaction_type = montecarlo_configuration.LINE_INTERACTION_TYPE

    if montecarlo_configuration.ENABLE_FULL_RELATIVITY:
        set_packet_props_full_relativity(r_packet, numba_model)
    else:
        set_packet_props_partial_relativity(r_packet, numba_model)
    r_packet.initialize_line_id(
        opacity_state,
        numba_model,
        montecarlo_configuration.ENABLE_FULL_RELATIVITY,
    )

    trace_vpacket_volley(
        r_packet,
        vpacket_collection,
        numba_radial_1d_geometry,
        numba_model,
        opacity_state,
        montecarlo_configuration.ENABLE_FULL_RELATIVITY,
        montecarlo_configuration.VPACKET_TAU_RUSSIAN,
        montecarlo_configuration.SURVIVAL_PROBABILITY,
        montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
    )

    if montecarlo_configuration.ENABLE_RPACKET_TRACKING:
        rpacket_tracker.track(r_packet)

    # this part of the code is temporary and will be better incorporated
    while r_packet.status == PacketStatus.IN_PROCESS:
        # Compute continuum quantities
        # trace packet (takes opacities)
        doppler_factor = get_doppler_factor(
            r_packet.r,
            r_packet.mu,
            numba_model.time_explosion,
            montecarlo_configuration.ENABLE_FULL_RELATIVITY,
        )

        comov_nu = r_packet.nu * doppler_factor
        chi_e = chi_electron_calculator(
            opacity_state, comov_nu, r_packet.current_shell_id
        )
        if montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED:
            (
                chi_bf_tot,
                chi_bf_contributions,
                current_continua,
                x_sect_bfs,
                chi_ff,
            ) = chi_continuum_calculator(
                opacity_state, comov_nu, r_packet.current_shell_id
            )
            chi_continuum = chi_e + chi_bf_tot + chi_ff

            escat_prob = chi_e / chi_continuum  # probability of e-scatter
            if montecarlo_configuration.ENABLE_FULL_RELATIVITY:
                chi_continuum *= doppler_factor
            distance, interaction_type, delta_shell = trace_packet(
                r_packet,
                numba_radial_1d_geometry,
                numba_model,
                opacity_state,
                estimators,
                chi_continuum,
                escat_prob,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.DISABLE_LINE_SCATTERING,
            )
            update_bound_free_estimators(
                comov_nu,
                r_packet.energy * doppler_factor,
                r_packet.current_shell_id,
                distance,
                estimators,
                opacity_state.t_electrons[r_packet.current_shell_id],
                x_sect_bfs,
                current_continua,
                opacity_state.bf_threshold_list_nu,
            )
        else:
            escat_prob = 1.0
            chi_continuum = chi_e
            if montecarlo_configuration.ENABLE_FULL_RELATIVITY:
                chi_continuum *= doppler_factor
            distance, interaction_type, delta_shell = trace_packet(
                r_packet,
                numba_radial_1d_geometry,
                numba_model,
                opacity_state,
                estimators,
                chi_continuum,
                escat_prob,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.DISABLE_LINE_SCATTERING,
            )

        # If continuum processes: update continuum estimators

        if interaction_type == InteractionType.BOUNDARY:
            move_r_packet(
                r_packet,
                distance,
                numba_model.time_explosion,
                estimators,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            move_packet_across_shell_boundary(
                r_packet, delta_shell, len(numba_radial_1d_geometry.r_inner)
            )

        elif interaction_type == InteractionType.LINE:
            r_packet.last_interaction_type = 2
            move_r_packet(
                r_packet,
                distance,
                numba_model.time_explosion,
                estimators,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            line_scatter(
                r_packet,
                numba_model.time_explosion,
                line_interaction_type,
                opacity_state,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
            )

        elif interaction_type == InteractionType.ESCATTERING:
            r_packet.last_interaction_type = 1

            move_r_packet(
                r_packet,
                distance,
                numba_model.time_explosion,
                estimators,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            thomson_scatter(
                r_packet,
                numba_model.time_explosion,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
            )
        elif (
            montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
            and interaction_type == InteractionType.CONTINUUM_PROCESS
        ):
            r_packet.last_interaction_type = InteractionType.CONTINUUM_PROCESS
            move_r_packet(
                r_packet,
                distance,
                numba_model.time_explosion,
                estimators,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            continuum_event(
                r_packet,
                numba_model.time_explosion,
                opacity_state,
                chi_bf_tot,
                chi_ff,
                chi_bf_contributions,
                current_continua,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                numba_model,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
                montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
            )
        else:
            pass
        if montecarlo_configuration.ENABLE_RPACKET_TRACKING:
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
        enable_full_relativity=False,
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
        enable_full_relativity=True,
    )

    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor
    r_packet.mu = (r_packet.mu + beta) / (1 + beta * r_packet.mu)
