from numba import njit

from tardis import constants as const
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacities import (
    chi_continuum_calculator,
    chi_electron_calculator,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
)
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bound_free,
)
from tardis.transport.montecarlo.interaction_event_callers import (
    continuum_event,
    line_scatter_event,
)
from tardis.transport.montecarlo.interaction_events import (
    thomson_scatter,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    move_packet_across_shell_boundary,
    move_r_packet,
    trace_packet,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
    RPacket,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@njit
def packet_propagation(
    r_packet: RPacket,
    numba_radial_1d_geometry: NumbaRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumba,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    estimators_continuum: EstimatorsContinuum,
    rpacket_tracker,  # Excluded from type hints as it might be different types
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """
    Execute Monte Carlo transport for a single radiative packet.

    This function performs the complete Monte Carlo transport simulation for a
    single r-packet, handling all interactions including line scattering,
    electron scattering, and continuum processes. The packet is traced through
    the ejecta until it escapes or is absorbed.

    Parameters
    ----------
    r_packet
        The radiative packet to transport through the ejecta.
    numba_radial_1d_geometry
        The spherically symmetric geometry of the supernova ejecta.
    time_explosion
        Time since explosion in seconds.
    opacity_state
        Current opacity state containing line and continuum opacities.
    estimators_bulk
        Monte Carlo estimators for cell-level bulk radiation field quantities.
    estimators_line
        Monte Carlo estimators for line-level radiation field quantities.
    estimators_continuum
        Monte Carlo estimators for continuum interaction quantities.
    vpacket_collection
        Collection for storing virtual packets when enabled.
    rpacket_tracker
        Tracker for recording packet interactions and trajectories.
    montecarlo_configuration
        Configuration parameters for the Monte Carlo simulation.

    Returns
    -------
    This function modifies the r_packet object in-place and updates
    estimators and collections. No return value.

    """
    line_interaction_type = montecarlo_configuration.LINE_INTERACTION_TYPE

    # IIP mode: always use full relativity
    set_packet_props_full_relativity(r_packet, time_explosion)
    r_packet.initialize_line_id(
        opacity_state,
        time_explosion,
        enable_full_relativity=True,
    )

    rpacket_tracker.track_boundary_event(
        r_packet, from_shell_id=-1, to_shell_id=0
    )

    # this part of the code is temporary and will be better incorporated
    while r_packet.status == PacketStatus.IN_PROCESS:
        # Compute continuum quantities
        # trace packet (takes opacities)
        doppler_factor = get_doppler_factor(
            r_packet.r,
            r_packet.mu,
            time_explosion,
            enable_full_relativity=True,
        )

        comov_nu = r_packet.nu * doppler_factor
        chi_e = chi_electron_calculator(
            opacity_state, comov_nu, r_packet.current_shell_id
        )
        # IIP mode: continuum processes always enabled
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
        # IIP mode: full relativity always enabled
        chi_continuum *= doppler_factor
        distance, interaction_type, delta_shell = trace_packet(
            r_packet,
            numba_radial_1d_geometry,
            time_explosion,
            opacity_state,
            estimators_line,
            chi_continuum,
            escat_prob,
            enable_full_relativity=True,
            disable_line_scattering=montecarlo_configuration.DISABLE_LINE_SCATTERING,
        )
        update_estimators_bound_free(
            comov_nu,
            r_packet.energy * doppler_factor,
            r_packet.current_shell_id,
            distance,
            estimators_continuum,
            opacity_state.t_electrons[r_packet.current_shell_id],
            x_sect_bfs,
            current_continua,
            opacity_state.bf_threshold_list_nu,
            chi_ff,
        )

        # Handle interaction types
        if interaction_type == InteractionType.BOUNDARY:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                enable_full_relativity=True,
            )
            rpacket_tracker.track_boundary_event(
                r_packet,
                r_packet.current_shell_id,
                r_packet.current_shell_id + delta_shell,
            )

            move_packet_across_shell_boundary(
                r_packet,
                delta_shell,
                len(numba_radial_1d_geometry.r_inner),
            )

        elif interaction_type == InteractionType.LINE:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                enable_full_relativity=True,
            )

            rpacket_tracker.track_line_interaction_before(r_packet)

            line_scatter_event(
                r_packet,
                time_explosion,
                line_interaction_type,
                opacity_state,
                enable_full_relativity=True,
            )
            rpacket_tracker.track_line_interaction_after(r_packet)

        elif interaction_type == InteractionType.ESCATTERING:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                enable_full_relativity=True,
            )
            rpacket_tracker.track_escattering_interaction_before(r_packet)
            thomson_scatter(
                r_packet,
                time_explosion,
                enable_full_relativity=True,
            )
            rpacket_tracker.track_escattering_interaction_after(r_packet)

        # IIP mode: continuum processes always enabled
        elif interaction_type == InteractionType.CONTINUUM_PROCESS:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                enable_full_relativity=True,
            )
            rpacket_tracker.track_continuum_interaction_before(r_packet)
            continuum_event(
                r_packet,
                time_explosion,
                opacity_state,
                chi_bf_tot,
                chi_ff,
                chi_bf_contributions,
                current_continua,
                enable_full_relativity=True,
            )

            rpacket_tracker.track_continuum_interaction_after(r_packet)

        else:
            # Handle any unrecognized interaction types
            rpacket_tracker.track_boundary_event(
                r_packet, from_shell_id=-1, to_shell_id=0
            )

    # Registering the final boundary interaction.
    # Only for RPacketTracker
    # This is required by the RPacketPlotter tool
    #
    # NOTE: This approach assumes packets are instantaneous and records
    # their final state immediately. In a future time-dependent implementation,
    # this may need to be modified to account for packet propagation time
    # and potentially delayed final state recording.
    # Track final packet state as simple boundary event (current shell -> current shell + 1)
    # This explicitly records the packet's exit from the simulation domain
    rpacket_tracker.track_boundary_event(
        r_packet,
        from_shell_id=r_packet.current_shell_id,
        to_shell_id=r_packet.current_shell_id + 1,
    )


@njit
def set_packet_props_partial_relativity(
    r_packet: RPacket, time_explosion: float
) -> None:
    """
    Set packet properties using partial relativistic corrections.

    This function applies inverse Doppler corrections to the packet frequency
    and energy based on partial relativistic treatment (first-order in v/c).

    Parameters
    ----------
    r_packet
        The radiative packet whose properties will be modified.
    time_explosion
        Time since explosion in seconds, used to calculate velocity.

    Returns
    -------
    Modifies r_packet.nu and r_packet.energy in-place.
    """
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        time_explosion,
        enable_full_relativity=False,
    )
    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor


@njit
def set_packet_props_full_relativity(
    r_packet: RPacket, time_explosion: float
) -> None:
    """
    Set packet properties using full relativistic corrections.

    This function applies inverse Doppler corrections to the packet frequency,
    energy, and direction cosine based on full relativistic treatment, including
    aberration effects on the direction cosine.

    Parameters
    ----------
    r_packet
        The radiative packet whose properties will be modified.
    time_explosion
        Time since explosion in seconds, used to calculate velocity.

    Returns
    -------
    Modifies r_packet.nu, r_packet.energy, and r_packet.mu in-place.
    """
    beta = (r_packet.r / time_explosion) / C_SPEED_OF_LIGHT

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r,
        r_packet.mu,
        time_explosion,
        enable_full_relativity=True,
    )

    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor
    r_packet.mu = (r_packet.mu + beta) / (1 + beta * r_packet.mu)
