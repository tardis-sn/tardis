"""Non-homologous mode rad packet transport - line-only without continuum processes."""

import numpy as np
from numba import njit

from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.frame_transformations import (
    get_doppler_factor_nonhomologous,
)
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_line_nonhomologous,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_bulk,
)
from tardis.transport.montecarlo.nonhomologous_grid import (
    piecewise_linear_dvdr,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
    RPacket,
)


@njit(**njit_dict_no_parallel)
def trace_packet(
    r_packet: RPacket,
    numba_radial_1d_geometry: NumbaNonhomologousRadial1DGeometry,
    opacity_state: OpacityStateNumba,
    estimators_line: EstimatorsLine,
    opacity_electron: float,
    enable_full_relativity: bool,
    disable_line_scattering: bool,
) -> tuple[float, InteractionType, int]:
    """
    Traces the RPacket through the ejecta and stops when an interaction happens.

    Non-homologous mode: only handles line interactions and electron scattering.

    Parameters
    ----------
    r_packet : RPacket
        The radiative packet being transported
    numba_radial_1d_geometry : NumbaNonhomologousRadial1DGeometry
        Radial 1D geometry of the model
    opacity_state : OpacityStateNumba
        Opacity state containing line list and tau sobolev
    estimators_line : EstimatorsLine
        Line-level radiation field estimators
    opacity_electron : float
        Electron scattering opacity
    enable_full_relativity : bool
        Flag to enable full relativistic calculations
    disable_line_scattering : bool
        Flag to disable line scattering

    Returns
    -------
    tuple[float, InteractionType, int]
        (distance, interaction_type, delta_shell)
    """
    r_inner = numba_radial_1d_geometry.r_inner[r_packet.current_shell_id]
    r_outer = numba_radial_1d_geometry.r_outer[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining taus
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    v, dvdr = piecewise_linear_dvdr(r_packet.r, r_packet.current_shell_id, numba_radial_1d_geometry)

    # defining start for line interaction
    # If redshifting, use next line and line list in order
    # If blueshifting, use previous line and reverse line list
    if dvdr >= 0.0:
        start_line_id = r_packet.next_line_id
        loop_lim, loop_direction = len(opacity_state.line_list_nu), 1
    else:
        start_line_id = r_packet.prev_line_id
        loop_lim, loop_direction = -1, -1

    distance_electron = tau_event / opacity_electron
    cur_line_id = start_line_id  # initializing varibale for Numba
    # - do not remove
    last_line_id = len(opacity_state.line_list_nu) - 1
    for cur_line_id in range(start_line_id, loop_lim, loop_direction):
        # Going through the lines
        nu_line = opacity_state.line_list_nu[cur_line_id]

        # Getting the tau for the next line
        tau_trace_line = opacity_state.tau_sobolev[cur_line_id, r_packet.current_shell_id]

        # Adding it to the tau_trace_line_combined
        tau_trace_line_combined += tau_trace_line

        # Calculating the distance until the current photons co-moving nu
        # redshifts to the line frequency
        is_last_line = cur_line_id == last_line_id

        distance_trace = calculate_distance_line_nonhomologous(
            r_packet,
            numba_radial_1d_geometry,
            nu_line,
        )

        # calculating the tau electron of how far the trace has progressed
        tau_trace_electron = opacity_electron * distance_trace

        # calculating the trace
        tau_trace_combined = tau_trace_line_combined + tau_trace_electron

        distance = min(distance_trace, distance_boundary, distance_electron)

        if distance_trace != 0:
            if (distance == distance_boundary) or (distance == distance_electron):
                if dvdr >= 0.0:
                    r_packet.next_line_id = cur_line_id
                    r_packet.prev_line_id = cur_line_id - 1
                else:
                    r_packet.next_line_id = cur_line_id + 1
                    r_packet.prev_line_id = cur_line_id
            if distance == distance_boundary:
                interaction_type = InteractionType.BOUNDARY  # BOUNDARY
                break
            if distance == distance_electron:
                interaction_type = InteractionType.ESCATTERING
                break

        # Updating the J_b_lu and E_dot_lu
        # This means we are still looking for line interaction and have not
        # been kicked out of the path by boundary or electron interaction

        # TODO:nonhomology - replace this function with a generalized version
        # Will need to pass geometry instead of t_explosion to get correct doppler factor
        #update_estimators_line(
        #    estimators_line,
        #    r_packet,
        #    cur_line_id,
        #    distance_trace,
        #    time_explosion,
        #    enable_full_relativity,
        #)
        # connor-mcclellan: here I reproduce the steps of the update_estimators_line -> calc_packet_energy call
        # stack, modifying the expressions to not use homologous expansion
        # This likely belongs somewhere else and can be moved in a future restructure, but for now
        # I'll keep in within the nonhomologous mode
        #
        # First step: getting the packet's new energy to use for the estimator update
        # Replaces the call to `calc_packet_energy` within `update_estimators_line`
        new_r = np.sqrt(
            r_packet.r * r_packet.r + distance_trace * distance_trace + 2.0 * r_packet.r * distance_trace * r_packet.mu
        )
        new_mu = (r_packet.mu * r_packet.r + distance_trace) / new_r
        new_v, _ = piecewise_linear_dvdr(new_r, r_packet.current_shell_id, numba_radial_1d_geometry)
        new_doppler_factor = (1.0 - new_v/C_SPEED_OF_LIGHT * new_mu)
        energy = r_packet.energy * new_doppler_factor

        # Second step: update the estimators
        # Replaces the call to `update_estimators_line`
        estimators_line.mean_intensity_blueward[
            cur_line_id, r_packet.current_shell_id
        ] += energy / r_packet.nu
        estimators_line.energy_deposition_line_rate[
            cur_line_id, r_packet.current_shell_id
        ] += energy

        if tau_trace_combined > tau_event and not disable_line_scattering:
            interaction_type = InteractionType.LINE  # Line
            if dvdr >= 0.0:
                r_packet.next_line_id = cur_line_id
                r_packet.prev_line_id = cur_line_id - 1
            else:
                r_packet.next_line_id = cur_line_id + 1
                r_packet.prev_line_id = cur_line_id
            distance = distance_trace
            break

        # Recalculating distance_electron using tau_event -
        # tau_trace_line_combined
        distance_electron = (tau_event - tau_trace_line_combined) / (
            opacity_electron
        )

    else:  # Executed when no break occurs in the for loop
        # We are beyond the line list now and the only next thing is to see
        # if we are interacting with the boundary or electron scattering
        if cur_line_id == (len(opacity_state.line_list_nu) - 1):
            # Treatment for last line
            cur_line_id += 1
        if distance_electron < distance_boundary:
            distance = distance_electron
            interaction_type = InteractionType.ESCATTERING
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    return distance, interaction_type, delta_shell


@njit(**njit_dict_no_parallel)
def move_r_packet(
    r_packet: RPacket,
    distance: float,
    geometry: NumbaNonhomologousRadial1DGeometry,
    estimators_bulk: EstimatorsBulk,
    enable_full_relativity: bool,
) -> None:
    """
    Move packet a distance and recalculate the new angle mu.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet object
    distance : float
        Distance to move in cm
    geometry : NumbaNonhomologousRadial1DGeometry
        Radius, velocity, and velocity gradient per cell.
    estimators_bulk : EstimatorsBulk
        Cell-level bulk radiation field estimators
    enable_full_relativity : bool
        Flag to enable full relativistic calculations
    """
    doppler_factor = get_doppler_factor_nonhomologous(
        r_packet.r,
        r_packet.mu,
        geometry,
        r_packet.current_shell_id,
        enable_full_relativity,
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
    packet : RPacket
        Radiative packet object
    delta_shell : int
        Change in shell index (+1 if moving outward or -1 if moving inward)
    no_of_shells : int
        Number of shells in TARDIS simulation
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id
