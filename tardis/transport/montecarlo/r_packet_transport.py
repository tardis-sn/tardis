import numpy as np
from numba import njit

import tardis.transport.montecarlo.configuration.montecarlo_globals as montecarlo_globals
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor,
)
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_line,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bulk,
    update_estimators_line,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
    RPacket,
)


@njit(**njit_dict_no_parallel)
def trace_packet(
    r_packet: RPacket,
    numba_radial_1d_geometry: NumbaRadial1DGeometry,
    time_explosion: float,
    opacity_state: OpacityStateNumba,
    estimators_line: EstimatorsLine,
    chi_continuum: float,
    escat_prob: float,
    enable_full_relativity: bool,
    disable_line_scattering: bool,
) -> tuple:
    """
    Traces the RPacket through the ejecta and stops when an interaction happens (heart of the calculation).

    Parameters
    ----------
    r_packet
        The radiative packet being transported
    numba_radial_1d_geometry
        Radial 1D geometry of the model
    time_explosion
        Time since explosion in seconds
    opacity_state
        Opacity state containing line list and tau sobolev
    estimators_line
        Line-level radiation field estimators
    chi_continuum
        Continuum opacity
    escat_prob
        Probability of electron scattering
    enable_full_relativity
        Flag to enable full relativistic calculations
    disable_line_scattering
        Flag to disable line scattering

    Returns
    -------
    tuple
        (distance, interaction_type, delta_shell)
    """
    r_inner = numba_radial_1d_geometry.r_inner[r_packet.current_shell_id]
    r_outer = numba_radial_1d_geometry.r_outer[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining start for line interaction
    start_line_id = r_packet.next_line_id

    # defining taus
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    # Calculating doppler factor
    doppler_factor = get_doppler_factor(
        r_packet.r,
        r_packet.mu,
        time_explosion,
        enable_full_relativity,
    )
    comov_nu = r_packet.nu * doppler_factor

    distance_continuum = tau_event / chi_continuum
    cur_line_id = start_line_id  # initializing varibale for Numba
    # - do not remove
    last_line_id = len(opacity_state.line_list_nu) - 1
    for cur_line_id in range(start_line_id, len(opacity_state.line_list_nu)):
        # Going through the lines
        nu_line = opacity_state.line_list_nu[cur_line_id]

        # Getting the tau for the next line
        tau_trace_line = opacity_state.tau_sobolev[
            cur_line_id, r_packet.current_shell_id
        ]

        # Adding it to the tau_trace_line_combined
        tau_trace_line_combined += tau_trace_line

        # Calculating the distance until the current photons co-moving nu
        # redshifts to the line frequency
        is_last_line = cur_line_id == last_line_id

        distance_trace = calculate_distance_line(
            r_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
            enable_full_relativity,
        )

        # calculating the tau continuum of how far the trace has progressed
        tau_trace_continuum = chi_continuum * distance_trace

        # calculating the trace
        tau_trace_combined = tau_trace_line_combined + tau_trace_continuum

        distance = min(distance_trace, distance_boundary, distance_continuum)

        if distance_trace != 0:
            if distance == distance_boundary:
                interaction_type = InteractionType.BOUNDARY  # BOUNDARY
                r_packet.next_line_id = cur_line_id
                break
            if distance == distance_continuum:
                if not montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
                    interaction_type = InteractionType.ESCATTERING
                else:
                    zrand = np.random.random()
                    if zrand < escat_prob:
                        interaction_type = InteractionType.ESCATTERING
                    else:
                        interaction_type = InteractionType.CONTINUUM_PROCESS
                r_packet.next_line_id = cur_line_id
                break

        # Updating the J_b_lu and E_dot_lu
        # This means we are still looking for line interaction and have not
        # been kicked out of the path by boundary or electron interaction

        update_estimators_line(
            estimators_line,
            r_packet,
            cur_line_id,
            distance_trace,
            time_explosion,
            enable_full_relativity,
        )

        if tau_trace_combined > tau_event and not disable_line_scattering:
            interaction_type = InteractionType.LINE  # Line
            r_packet.next_line_id = cur_line_id
            distance = distance_trace
            break

        # Recalculating distance_continuum using tau_event -
        # tau_trace_line_combined
        # I don't think this needs to be updated
        # since tau_event is already the result of the integral
        # from the initial line
        distance_continuum = (tau_event - tau_trace_line_combined) / (
            chi_continuum
        )

    else:  # Executed when no break occurs in the for loop
        # We are beyond the line list now and the only next thing is to see
        # if we are interacting with the boundary or electron scattering
        if cur_line_id == (len(opacity_state.line_list_nu) - 1):
            # Treatment for last line
            cur_line_id += 1
        if distance_continuum < distance_boundary:
            distance = distance_continuum
            if not montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
                interaction_type = InteractionType.ESCATTERING
            else:
                zrand = np.random.random()
                if zrand < escat_prob:
                    interaction_type = InteractionType.ESCATTERING
                else:
                    interaction_type = InteractionType.CONTINUUM_PROCESS
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    return distance, interaction_type, delta_shell


@njit(**njit_dict_no_parallel)
def move_r_packet(
    r_packet: RPacket,
    distance: float,
    time_explosion: float,
    estimators_bulk: EstimatorsBulk,
    enable_full_relativity: bool,
) -> None:
    """
    Move packet a distance and recalculate the new angle mu.

    Parameters
    ----------
    r_packet
        Radiative packet object
    distance
        Distance to move in cm
    time_explosion
        Time since explosion in seconds
    estimators_bulk
        Cell-level bulk radiation field estimators
    enable_full_relativity
        Flag to enable full relativistic calculations
    """
    doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion, enable_full_relativity
    )

    r = r_packet.r
    if distance > 0.0:
        new_r = np.sqrt(
            r * r + distance * distance + 2.0 * r * distance * r_packet.mu
        )
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

        comov_nu = r_packet.nu * doppler_factor
        comov_energy = r_packet.energy * doppler_factor

        # Account for length contraction
        if enable_full_relativity:
            distance *= doppler_factor

        update_estimators_bulk(
            r_packet, distance, estimators_bulk, comov_nu, comov_energy
        )


@njit(**njit_dict_no_parallel)
def move_packet_across_shell_boundary(
    packet: RPacket, delta_shell: int, no_of_shells: int
) -> None:
    """
    Move packet across shell boundary - realizing if we are still in the simulation or have
    moved out through the inner boundary or outer boundary and updating packet status.

    Parameters
    ----------
    packet
        Radiative packet object
    delta_shell
        Change in shell index (+1 if moving outward or -1 if moving inward)
    no_of_shells
        Number of shells in TARDIS simulation
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id
