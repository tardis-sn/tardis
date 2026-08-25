"""Classic mode packet transport without continuum processes."""

from numba import njit

from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacities import chi_electron_calculator
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import get_doppler_factor
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.modes.homologous_packet_propagation import (
    handle_common_interaction,
    packet_propagation_homologous,
    set_packet_props_full_relativity,
    set_packet_props_partial_relativity,
)
from tardis.transport.montecarlo.modes.homologous_rad_packet_transport import (
    trace_packet,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)
from tardis.transport.montecarlo.packets.virtual_packet import (
    trace_vpacket_volley,
)


@njit
def packet_propagation(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumba,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    vpacket_collection: VPacketCollection,
    rpacket_tracker: TrackerFull | TrackerLastInteraction,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """Execute line-only Monte Carlo transport for one radiative packet."""
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
        geometry,
        time_explosion,
        opacity_state,
        montecarlo_configuration.ENABLE_FULL_RELATIVITY,
        montecarlo_configuration.VPACKET_TAU_RUSSIAN,
        montecarlo_configuration.SURVIVAL_PROBABILITY,
    )

    packet_propagation_homologous(
        r_packet,
        geometry,
        time_explosion,
        opacity_state,
        estimators_bulk,
        estimators_line,
        vpacket_collection,
        rpacket_tracker,
        montecarlo_configuration,
        packet_propagation_step,
    )


@njit
def packet_propagation_step(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumba,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    vpacket_collection: VPacketCollection,
    rpacket_tracker: TrackerFull | TrackerLastInteraction,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """Propagate a classic packet through one interaction."""
    enable_full_relativity = montecarlo_configuration.ENABLE_FULL_RELATIVITY
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

    interaction_handled = handle_common_interaction(
        r_packet,
        geometry,
        time_explosion,
        opacity_state,
        estimators_bulk,
        rpacket_tracker,
        distance,
        interaction_type,
        delta_shell,
        montecarlo_configuration.LINE_INTERACTION_TYPE,
        enable_full_relativity,
    )
    if not interaction_handled:
        rpacket_tracker.track_boundary_event(
            r_packet, from_shell_id=-1, to_shell_id=0
        )
    elif interaction_type != InteractionType.BOUNDARY:
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
