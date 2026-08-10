"""Non-homologous mode rad packet transport - line-only without continuum processes."""

import numpy as np
from numba import njit

from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_line_nonhomologous,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    MISS_DISTANCE,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.utils import MonteCarloException


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
    v_inner = numba_radial_1d_geometry.v_inner[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining taus
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    dvdr = numba_radial_1d_geometry.velocity_gradient[
        r_packet.current_shell_id
    ]
    distance_electron = tau_event / opacity_electron
    distance_trace_previous = 0.0
    while True:
        distance_trace = MISS_DISTANCE
        cur_line_id = -1
        for line_id in range(len(opacity_state.line_list_nu)):
            nu_line = opacity_state.line_list_nu[line_id]
            distance_trace_line = calculate_distance_line_nonhomologous(
                r_packet,
                numba_radial_1d_geometry,
                nu_line,
                distance_trace_previous,
            )
            if distance_trace_line < distance_trace:
                distance_trace = distance_trace_line
                cur_line_id = line_id

        if distance_electron < min(distance_trace, distance_boundary):
            distance = distance_electron
            interaction_type = InteractionType.ESCATTERING
            break
        if distance_boundary <= distance_trace:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY
            break

        # Updating the J_b_lu and E_dot_lu
        # This means we are still looking for line interaction and have not
        # been kicked out of the path by boundary or electron interaction

        # connor-mcclellan: some hardcoded overrides here until update_estimators_line
        # is updated and generalized to support nonhomologous geometry

        # Get the packet's new energy to use for the estimator update
        # Replaces the call to `calc_packet_energy` within `update_estimators_line`
        new_r = np.sqrt(
            r_packet.r * r_packet.r
            + distance_trace * distance_trace
            + 2.0 * r_packet.r * distance_trace * r_packet.mu
        )
        new_mu = (r_packet.mu * r_packet.r + distance_trace) / new_r
        dvdr = numba_radial_1d_geometry.velocity_gradient[
            r_packet.current_shell_id
        ]
        new_v = v_inner + dvdr * (new_r - r_inner)
        new_doppler_factor = 1.0 - new_v / C_SPEED_OF_LIGHT * new_mu
        energy = r_packet.energy * new_doppler_factor
        projected_velocity_gradient = (
            new_mu * new_mu * dvdr
            + (1.0 - new_mu * new_mu) * new_v / new_r
        )
        if projected_velocity_gradient == 0.0:
            raise MonteCarloException(
                "Sobolev optical depth is singular at the line resonance."
            )

        tau_trace_line = opacity_state.sobolev_line_strength[
            cur_line_id, r_packet.current_shell_id
        ]
        if tau_trace_line == 0.0:
            tau_trace_line = (
                opacity_state.tau_sobolev[
                    cur_line_id, r_packet.current_shell_id
                ]
                * abs(dvdr)
            )
        if disable_line_scattering:
            tau_trace_line = 0.0
        else:
            tau_trace_line /= abs(projected_velocity_gradient)

        tau_trace_electron = opacity_electron * distance_trace
        tau_trace_combined = (
            tau_trace_line_combined
            + tau_trace_line
            + tau_trace_electron
        )

        # Update the estimators
        # Replaces the call to `update_estimators_line`
        estimators_line.mean_intensity_blueward[
            cur_line_id, r_packet.current_shell_id
        ] += energy / r_packet.nu / abs(projected_velocity_gradient)
        estimators_line.energy_deposition_line_rate[
            cur_line_id, r_packet.current_shell_id
        ] += energy

        if tau_trace_combined > tau_event and not disable_line_scattering:
            interaction_type = InteractionType.LINE  # Line
            r_packet.next_line_id = cur_line_id
            r_packet.prev_line_id = cur_line_id - 1
            distance = distance_trace
            break

        tau_trace_line_combined += tau_trace_line
        distance_electron = (tau_event - tau_trace_line_combined) / (
            opacity_electron
        )
        distance_trace_previous = distance_trace

    return distance, interaction_type, delta_shell
