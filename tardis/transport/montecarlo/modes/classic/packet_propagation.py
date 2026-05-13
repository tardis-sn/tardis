"""Classic mode packet transport - line-only without continuum processes."""

from numba import njit

from tardis import constants as const
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacities import chi_electron_calculator
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
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.interaction_event_callers import (
    line_scatter_event,
)
from tardis.transport.montecarlo.interaction_events import (
    thomson_scatter,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
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
from tardis.transport.montecarlo.packets.virtual_packet import (
    trace_vpacket_volley,
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
    vpacket_collection: VPacketCollection,
    rpacket_tracker,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """
    Execute Monte Carlo transport for a single radiative packet in classic mode.

    Classic mode performs line-only transport without continuum processes.
    Only electron scattering and line interactions are handled.

    Parameters
    ----------
    r_packet : RPacket
        The radiative packet to transport through the ejecta.
    numba_radial_1d_geometry : NumbaRadial1DGeometry
        The spherically symmetric geometry of the supernova ejecta.
    time_explosion : float
        Time since explosion in seconds.
    opacity_state : OpacityStateNumba
        Current opacity state containing line opacities.
    estimators_bulk : EstimatorsBulk
        Monte Carlo estimators for cell-level bulk radiation field quantities.
    estimators_line : EstimatorsLine
        Monte Carlo estimators for line-level radiation field quantities.
    vpacket_collection : VPacketCollection
        Collection for storing virtual packets when enabled.
    rpacket_tracker
        Tracker for recording packet interactions and trajectories.
    montecarlo_configuration : MonteCarloConfiguration
        Configuration parameters for the Monte Carlo simulation.

    Returns
    -------
    None
        This function modifies the r_packet object in-place and updates
        estimators and collections.
    """
    line_interaction_type = montecarlo_configuration.LINE_INTERACTION_TYPE

    if montecarlo_configuration.ENABLE_FULL_RELATIVITY:
        set_packet_props_full_relativity(r_packet, time_explosion)
    else:
        set_packet_props_partial_relativity(r_packet, time_explosion)
    r_packet.initialize_line_id(
        opacity_state,
        time_explosion,
        montecarlo_configuration.ENABLE_FULL_RELATIVITY,
    )

    trace_vpacket_volley(
        r_packet,
        vpacket_collection,
        numba_radial_1d_geometry,
        time_explosion,
        opacity_state,
        montecarlo_configuration.ENABLE_FULL_RELATIVITY,
        montecarlo_configuration.VPACKET_TAU_RUSSIAN,
        montecarlo_configuration.SURVIVAL_PROBABILITY,
    )

    rpacket_tracker.track_boundary_event(
        r_packet, from_shell_id=-1, to_shell_id=0
    )

    # this part of the code is temporary and will be better incorporated
    while r_packet.status == PacketStatus.IN_PROCESS:
        # Compute electron scattering opacity
        doppler_factor = get_doppler_factor(
            r_packet.r,
            r_packet.mu,
            time_explosion,
            montecarlo_configuration.ENABLE_FULL_RELATIVITY,
        )

        comov_nu = r_packet.nu * doppler_factor
        opacity_electron = chi_electron_calculator(
            opacity_state, comov_nu, r_packet.current_shell_id
        )

        if montecarlo_configuration.ENABLE_FULL_RELATIVITY:
            opacity_electron *= doppler_factor

        distance, interaction_type, delta_shell = trace_packet(
            r_packet,
            numba_radial_1d_geometry,
            time_explosion,
            opacity_state,
            estimators_line,
            opacity_electron,
            montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            montecarlo_configuration.DISABLE_LINE_SCATTERING,
        )

        if interaction_type == InteractionType.BOUNDARY:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
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
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )

            rpacket_tracker.track_line_interaction_before(r_packet)

            line_scatter_event(
                r_packet,
                time_explosion,
                line_interaction_type,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            rpacket_tracker.track_line_interaction_after(r_packet)
            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                time_explosion,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
            )

        elif interaction_type == InteractionType.ESCATTERING:
            move_r_packet(
                r_packet,
                distance,
                time_explosion,
                estimators_bulk,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            rpacket_tracker.track_escattering_interaction_before(r_packet)
            thomson_scatter(
                r_packet,
                time_explosion,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
            )
            rpacket_tracker.track_escattering_interaction_after(r_packet)

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                numba_radial_1d_geometry,
                time_explosion,
                opacity_state,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
            )
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
