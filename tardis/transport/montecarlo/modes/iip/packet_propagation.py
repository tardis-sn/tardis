"""IIP packet transport with continuum processes."""

from numba import njit

from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacities import (
    chi_continuum_calculator,
    chi_electron_calculator,
)
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.transport.frame_transformations import get_doppler_factor
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
)
from tardis.transport.montecarlo.modes.homologous_packet_propagation import (
    handle_common_interaction,
    packet_propagation_homologous,
    set_packet_props_full_relativity,
)
from tardis.transport.montecarlo.modes.homologous_rad_packet_transport import (
    trace_packet,
)
from tardis.transport.montecarlo.packets.movement import move_r_packet
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


@njit
def packet_propagation(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumbaIIP,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    estimators_continuum: EstimatorsContinuum,
    rpacket_tracker: TrackerFull | TrackerLastInteraction,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """Execute continuum Monte Carlo transport for one radiative packet."""
    set_packet_props_full_relativity(r_packet, time_explosion)
    r_packet.initialize_line_id(
        opacity_state,
        time_explosion,
        enable_full_relativity=True,
    )

    packet_propagation_homologous(
        r_packet,
        geometry,
        time_explosion,
        opacity_state,
        estimators_bulk,
        estimators_line,
        estimators_continuum,
        rpacket_tracker,
        montecarlo_configuration,
        packet_propagation_step,
    )


@njit
def packet_propagation_step(
    r_packet: RPacket,
    geometry: NumbaHomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumbaIIP,
    estimators_bulk: EstimatorsBulk,
    estimators_line: EstimatorsLine,
    estimators_continuum: EstimatorsContinuum,
    rpacket_tracker: TrackerFull | TrackerLastInteraction,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """Propagate an IIP packet through one interaction."""
    velocity = geometry.get_velocity(r_packet.r, r_packet.current_shell_id)
    doppler_factor = get_doppler_factor(
        velocity,
        r_packet.mu,
        enable_full_relativity=True,
    )
    comov_nu = r_packet.nu * doppler_factor
    chi_e = chi_electron_calculator(
        opacity_state, comov_nu, r_packet.current_shell_id
    )
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
    escat_prob = chi_e / chi_continuum
    chi_continuum *= doppler_factor

    distance, interaction_type, delta_shell = trace_packet(
        r_packet,
        geometry,
        time_explosion,
        opacity_state,
        estimators_line,
        chi_continuum,
        escat_prob,
        True,
        enable_full_relativity=True,
        disable_line_scattering=(
            montecarlo_configuration.DISABLE_LINE_SCATTERING
        ),
    )
    update_estimators_bound_free(
        comov_nu,
        r_packet.energy * doppler_factor,
        r_packet.current_shell_id,
        distance * doppler_factor,
        estimators_continuum,
        opacity_state.t_electrons[r_packet.current_shell_id],
        x_sect_bfs,
        current_continua,
        opacity_state.bf_threshold_list_nu,
        chi_ff * doppler_factor,
    )

    if handle_common_interaction(
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
        True,
    ):
        return

    if interaction_type == InteractionType.CONTINUUM_PROCESS:
        move_r_packet(
            r_packet,
            distance,
            geometry,
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
        rpacket_tracker.track_boundary_event(
            r_packet, from_shell_id=-1, to_shell_id=0
        )
