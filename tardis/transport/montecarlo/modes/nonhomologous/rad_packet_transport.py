"""Non-homologous mode rad packet transport - line-only without continuum processes."""

import numpy as np
import numpy.typing as npt
from numba import njit

from tardis.model.geometry.radial1d import (
    NumbaRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.geometry.calculate_distances import (
    calculate_comoving_frequency_nonhomologous,
    calculate_distance_boundary,
    calculate_distance_line,
    calculate_projected_gradient_zero_distances,
    get_line_id_range_nonhomologous,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.utils import MonteCarloException


@njit(**njit_dict_no_parallel)
def update_line_ids(
    r_packet: RPacket,
    line_list_nu: npt.NDArray[np.float64],
    comov_nu: float,
) -> None:
    """Update the line IDs bracketing a packet's comoving frequency.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet whose line IDs are updated.
    line_list_nu : numpy.ndarray
        Line frequencies in descending order.
    comov_nu : float
        Packet frequency in the comoving frame.
    """
    next_line_id = len(line_list_nu) - np.searchsorted(
        line_list_nu[::-1], comov_nu
    )
    r_packet.next_line_id = next_line_id
    r_packet.prev_line_id = next_line_id - 1


@njit(**njit_dict_no_parallel)
def trace_packet(
    r_packet: RPacket,
    numba_radial_1d_geometry: NumbaRadial1DGeometry,
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
    numba_radial_1d_geometry : NumbaRadial1DGeometry
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
    (
        first_gradient_zero_distance,
        second_gradient_zero_distance,
        gradient_zero_count,
    ) = calculate_projected_gradient_zero_distances(
        r_packet, numba_radial_1d_geometry, distance_boundary
    )

    interval_start = 0.0
    for interval_id in range(gradient_zero_count + 1):
        if interval_id == 0 and gradient_zero_count > 0:
            interval_end = first_gradient_zero_distance
        elif interval_id == 1 and gradient_zero_count > 1:
            interval_end = second_gradient_zero_distance
        else:
            interval_end = distance_boundary

        comov_nu_start = calculate_comoving_frequency_nonhomologous(
            r_packet, numba_radial_1d_geometry, interval_start
        )
        comov_nu_end = calculate_comoving_frequency_nonhomologous(
            r_packet, numba_radial_1d_geometry, interval_end
        )
        (
            start_line_id,
            stop_line_id,
            line_id_step,
        ) = get_line_id_range_nonhomologous(
            opacity_state.line_list_nu,
            comov_nu_start,
            comov_nu_end,
        )

        for cur_line_id in range(
            start_line_id, stop_line_id, line_id_step
        ):
            distance_trace = calculate_distance_line(
                r_packet,
                numba_radial_1d_geometry,
                opacity_state.line_list_nu[cur_line_id],
                interval_start,
                interval_end,
            )
            if distance_trace > interval_end:
                continue

            if distance_electron < distance_trace:
                comov_nu_event = calculate_comoving_frequency_nonhomologous(
                    r_packet,
                    numba_radial_1d_geometry,
                    distance_electron,
                )
                update_line_ids(
                    r_packet,
                    opacity_state.line_list_nu,
                    comov_nu_event,
                )
                return (
                    distance_electron,
                    InteractionType.ESCATTERING,
                    delta_shell,
                )
            if distance_boundary <= distance_trace:
                update_line_ids(
                    r_packet,
                    opacity_state.line_list_nu,
                    comov_nu_end,
                )
                return (
                    distance_boundary,
                    InteractionType.BOUNDARY,
                    delta_shell,
                )

            # connor-mcclellan: some hardcoded overrides here until
            # update_estimators_line is generalized for nonhomologous geometry
            new_r = np.sqrt(
                r_packet.r * r_packet.r
                + distance_trace * distance_trace
                + 2.0 * r_packet.r * distance_trace * r_packet.mu
            )
            new_mu = (r_packet.mu * r_packet.r + distance_trace) / new_r
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

            tau_trace_line = (
                opacity_state.sobolev_optical_depth_coefficient[
                    cur_line_id, r_packet.current_shell_id
                ]
            )
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

            estimators_line.mean_intensity_blueward[
                cur_line_id, r_packet.current_shell_id
            ] += energy / r_packet.nu / abs(projected_velocity_gradient)
            estimators_line.energy_deposition_line_rate[
                cur_line_id, r_packet.current_shell_id
            ] += energy

            if (
                tau_trace_combined > tau_event
                and not disable_line_scattering
            ):
                r_packet.next_line_id = cur_line_id
                r_packet.prev_line_id = cur_line_id - 1
                return distance_trace, InteractionType.LINE, delta_shell

            tau_trace_line_combined += tau_trace_line
            distance_electron = (
                tau_event - tau_trace_line_combined
            ) / opacity_electron

        if distance_electron < interval_end:
            comov_nu_event = calculate_comoving_frequency_nonhomologous(
                r_packet,
                numba_radial_1d_geometry,
                distance_electron,
            )
            update_line_ids(
                r_packet,
                opacity_state.line_list_nu,
                comov_nu_event,
            )
            return (
                distance_electron,
                InteractionType.ESCATTERING,
                delta_shell,
            )
        interval_start = interval_end

    update_line_ids(
        r_packet,
        opacity_state.line_list_nu,
        comov_nu_end,
    )
    return distance_boundary, InteractionType.BOUNDARY, delta_shell
