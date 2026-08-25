"""Shared packet propagation for homologous transport modes."""

from collections.abc import Callable

from numba import njit

from tardis import constants as const
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.transport.frame_transformations import (
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
from tardis.transport.montecarlo.interaction_event_callers import (
    line_scatter_event,
)
from tardis.transport.montecarlo.interaction_events import thomson_scatter
from tardis.transport.montecarlo.packets.movement import (
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
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)

HomologousOpacityState = OpacityStateNumba | OpacityStateNumbaIIP
HomologousModeState = VPacketCollection | EstimatorsContinuum
RPacketTracker = TrackerFull | TrackerLastInteraction

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@njit
def packet_propagation_homologous(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: HomologousOpacityState,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    mode_state: HomologousModeState,
    rpacket_tracker: RPacketTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    packet_propagation_step: Callable[..., None],
) -> None:
    """Propagate a packet using a compiled mode-specific transport step."""
    rpacket_tracker.track_boundary_event(
        r_packet, from_shell_id=-1, to_shell_id=0
    )

    while r_packet.status == PacketStatus.IN_PROCESS:
        packet_propagation_step(
            r_packet,
            geometry,
            time_explosion,
            opacity_state,
            estimators_bulk,
            estimators_line,
            mode_state,
            rpacket_tracker,
            montecarlo_configuration,
        )

    rpacket_tracker.track_boundary_event(
        r_packet,
        from_shell_id=r_packet.current_shell_id,
        to_shell_id=r_packet.current_shell_id + 1,
    )


@njit
def handle_common_interaction(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: HomologousOpacityState,
    estimators_bulk: EstimatorsBulk,
    rpacket_tracker: RPacketTracker,
    distance: float,
    interaction_type: int,
    delta_shell: int,
    line_interaction_type: int,
    enable_full_relativity: bool,
) -> bool:
    """Handle an interaction shared by classic and IIP transport."""
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
    else:
        return False

    return True


@njit
def set_packet_props_partial_relativity(
    r_packet: RPacket, time_explosion: float
) -> None:
    """Transform initial packet frequency and energy to the lab frame."""
    velocity = r_packet.r / time_explosion
    inverse_doppler_factor = get_inverse_doppler_factor(
        velocity,
        r_packet.mu,
        enable_full_relativity=False,
    )
    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor


@njit
def set_packet_props_full_relativity(
    r_packet: RPacket, time_explosion: float
) -> None:
    """Transform initial packet properties to the lab frame relativistically."""
    velocity = r_packet.r / time_explosion
    beta = velocity / C_SPEED_OF_LIGHT
    inverse_doppler_factor = get_inverse_doppler_factor(
        velocity,
        r_packet.mu,
        enable_full_relativity=True,
    )

    r_packet.nu *= inverse_doppler_factor
    r_packet.energy *= inverse_doppler_factor
    r_packet.mu = (r_packet.mu + beta) / (1 + beta * r_packet.mu)
