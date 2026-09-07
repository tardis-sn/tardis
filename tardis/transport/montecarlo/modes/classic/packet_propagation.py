"""Classic mode packet transport - line-only without continuum processes."""

from numba import njit

from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacities import chi_electron_calculator
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor,
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
from tardis.transport.montecarlo.modes.homologous_rad_packet_transport import (
    trace_packet,
)
from tardis.transport.montecarlo.packets.movement import (
    initialize_packet_frame,
    move_packet_across_shell_boundary,
    move_r_packet,
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


@njit
def packet_propagation(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    opacity_state: OpacityStateNumba,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    vpacket_collection: VPacketCollection,
    rpacket_tracker: object,
    montecarlo_configuration: MonteCarloConfiguration,
    enable_full_relativity: bool,
) -> None:
    """
    Execute Monte Carlo transport for a single radiative packet in classic mode.

    Classic mode performs line-only transport without continuum processes.
    Only electron scattering and line interactions are handled.

    Parameters
    ----------
    r_packet : RPacket
        The radiative packet to transport through the ejecta.
    geometry : NumbaHomologousRadial1DGeometry
        The spherically symmetric geometry of the supernova ejecta.
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
    enable_full_relativity : bool
        Whether to use TARDIS's existing full-relativity branch.

    """
    time_explosion = geometry.time_explosion
    line_interaction_type = montecarlo_configuration.LINE_INTERACTION_TYPE

    initialize_packet_frame(
        r_packet,
        geometry,
        enable_full_relativity,
    )
    r_packet.initialize_line_id(
        opacity_state,
        time_explosion,
        enable_full_relativity,
    )

    trace_vpacket_volley(
        r_packet,
        vpacket_collection,
        geometry,
        time_explosion,
        opacity_state,
        enable_full_relativity,
        montecarlo_configuration.VPACKET_TAU_RUSSIAN,
        montecarlo_configuration.SURVIVAL_PROBABILITY,
    )

    rpacket_tracker.track_boundary_event(
        r_packet, from_shell_id=-1, to_shell_id=0
    )

    # this part of the code is temporary and will be better incorporated
    while r_packet.status == PacketStatus.IN_PROCESS:
        # Compute electron scattering opacity
        velocity = geometry.get_velocity(r_packet.r, r_packet.current_shell_id)
        doppler_factor = get_doppler_factor(
            velocity,
            r_packet.mu,
            enable_full_relativity,
        )

        comov_nu = r_packet.nu * doppler_factor
        opacity_electron = chi_electron_calculator(
            opacity_state, comov_nu, r_packet.current_shell_id
        )

        if enable_full_relativity:
            opacity_electron *= doppler_factor

        distance, interaction_type, delta_shell = trace_packet(
            r_packet,
            geometry,
            time_explosion,
            opacity_state,
            estimators_line,
            opacity_electron,
            1.0,
            False,
            enable_full_relativity,
            montecarlo_configuration.DISABLE_LINE_SCATTERING,
        )

        if interaction_type == InteractionType.BOUNDARY:
            move_r_packet(
                r_packet,
                distance,
                geometry,
                estimators_bulk,
                enable_full_relativity,
            )
            rpacket_tracker.track_boundary_event(
                r_packet,
                r_packet.current_shell_id,
                r_packet.current_shell_id + delta_shell,
            )

            move_packet_across_shell_boundary(
                r_packet,
                delta_shell,
                len(geometry.r_inner),
            )

        elif interaction_type == InteractionType.LINE:
            move_r_packet(
                r_packet,
                distance,
                geometry,
                estimators_bulk,
                enable_full_relativity,
            )

            rpacket_tracker.track_line_interaction_before(r_packet)

            line_scatter_event(
                r_packet,
                time_explosion,
                line_interaction_type,
                opacity_state,
                enable_full_relativity,
            )
            rpacket_tracker.track_line_interaction_after(r_packet)
            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                geometry,
                time_explosion,
                opacity_state,
                enable_full_relativity,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
            )

        elif interaction_type == InteractionType.ESCATTERING:
            move_r_packet(
                r_packet,
                distance,
                geometry,
                estimators_bulk,
                enable_full_relativity,
            )
            rpacket_tracker.track_escattering_interaction_before(r_packet)
            thomson_scatter(
                r_packet,
                time_explosion,
                enable_full_relativity,
            )
            rpacket_tracker.track_escattering_interaction_after(r_packet)

            trace_vpacket_volley(
                r_packet,
                vpacket_collection,
                geometry,
                time_explosion,
                opacity_state,
                enable_full_relativity,
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
