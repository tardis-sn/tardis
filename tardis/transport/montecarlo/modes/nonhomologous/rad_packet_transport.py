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
    calculate_distance_line_nonhomologous,
    calculate_packet_position,
    calculate_projected_gradient_zero_distances,
    calculate_projected_velocity_gradient,
    get_line_id_range_nonhomologous,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    CLOSE_LINE_THRESHOLD,
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
def update_line_ids(
    r_packet: RPacket,
    line_list_nu: npt.NDArray[np.float64],
    comov_nu: float,
    active_line_id_step: int = 0,
) -> None:
    """Initialize eligible line IDs or refresh the opposite search cursor.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet whose line IDs are updated.
    line_list_nu : numpy.ndarray
        Line frequencies in descending order.
    comov_nu : float
        Packet frequency in the comoving frame.
    active_line_id_step : int, optional
        Preserve the redward cursor for ``+1`` or the blueward cursor for
        ``-1``. Zero initializes both cursors for a new frequency history.
    """
    frequency_padding = CLOSE_LINE_THRESHOLD * abs(r_packet.nu)
    ascending_line_nu = line_list_nu[::-1]
    line_count = len(line_list_nu)
    if active_line_id_step <= 0:
        r_packet.next_line_id = line_count - np.searchsorted(
            ascending_line_nu, comov_nu + frequency_padding, side="right"
        )
    if active_line_id_step >= 0:
        r_packet.prev_line_id = line_count - np.searchsorted(
            ascending_line_nu, comov_nu - frequency_padding, side="left"
        ) - 1


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
        Distance, interaction type, and shell-index change. The packet's
        directional line cursors are updated in place.
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

    dvdr = numba_radial_1d_geometry.velocity_gradient[r_packet.current_shell_id]
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
            CLOSE_LINE_THRESHOLD * abs(r_packet.nu) if interval_id == 0 else 0.0,
        )
        # Only the first interval inherits exclusions from the preceding trace.
        # A later reversal makes previously crossed lines eligible again.
        if interval_id == 0:
            if line_id_step == 1:
                start_line_id = max(start_line_id, r_packet.next_line_id)
            else:
                start_line_id = min(start_line_id, r_packet.prev_line_id)
        else:
            if line_id_step == 1:
                r_packet.next_line_id = start_line_id
            else:
                r_packet.prev_line_id = start_line_id

        for cur_line_id in range(start_line_id, stop_line_id, line_id_step):
            line_frequency = opacity_state.line_list_nu[cur_line_id]
            if interval_id == 0 and abs(
                (comov_nu_start - line_frequency) / r_packet.nu
            ) < CLOSE_LINE_THRESHOLD:
                # Eligible close lines resonate at the starting point;
                # processed line IDs have already been excluded above.
                distance_trace = 0.0
            else:
                distance_trace = calculate_distance_line_nonhomologous(
                    r_packet,
                    numba_radial_1d_geometry,
                    line_frequency,
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
                if distance_electron > 0.0:
                    update_line_ids(
                        r_packet,
                        opacity_state.line_list_nu,
                        comov_nu_event,
                        line_id_step,
                    )
                return (
                    distance_electron,
                    InteractionType.ESCATTERING,
                    delta_shell,
                )
            if distance_boundary <= distance_trace:
                if distance_boundary > 0.0:
                    update_line_ids(
                        r_packet,
                        opacity_state.line_list_nu,
                        comov_nu_end,
                        line_id_step,
                    )
                return (
                    distance_boundary,
                    InteractionType.BOUNDARY,
                    delta_shell,
                )

            # Line interaction
            (
                resonance_radius,
                resonance_direction_cosine,
            ) = calculate_packet_position(
                r_packet.r, r_packet.mu, distance_trace
            )
            velocity_at_resonance = numba_radial_1d_geometry.get_velocity(
                resonance_radius, r_packet.current_shell_id
            )
            projected_velocity_gradient = calculate_projected_velocity_gradient(
                resonance_radius,
                resonance_direction_cosine,
                velocity_at_resonance,
                dvdr,
            )
            new_doppler_factor = (
                1.0
                - velocity_at_resonance
                / C_SPEED_OF_LIGHT
                * resonance_direction_cosine
            )
            energy = r_packet.energy * new_doppler_factor
            if projected_velocity_gradient == 0.0:
                raise MonteCarloException(
                    "Sobolev optical depth is singular at the line resonance."
                )

            tau_trace_line = opacity_state.sobolev_optical_depth_coefficient[
                cur_line_id, r_packet.current_shell_id
            ]
            if tau_trace_line == 0.0:
                tau_trace_line = opacity_state.tau_sobolev[
                    cur_line_id, r_packet.current_shell_id
                ] * abs(dvdr)
            if disable_line_scattering:
                tau_trace_line = 0.0
            else:
                tau_trace_line /= abs(projected_velocity_gradient)

            tau_trace_electron = opacity_electron * distance_trace
            tau_trace_combined = (
                tau_trace_line_combined + tau_trace_line + tau_trace_electron
            )

            # Nonhomologous transport applies the local projected gradient at
            # each resonance; homologous transport can normalize globally.
            estimators_line.mean_intensity_blueward[
                cur_line_id, r_packet.current_shell_id
            ] += energy / r_packet.nu / abs(projected_velocity_gradient)
            estimators_line.energy_deposition_line_rate[
                cur_line_id, r_packet.current_shell_id
            ] += energy
            if line_id_step == 1:
                r_packet.next_line_id = cur_line_id + 1
            else:
                r_packet.prev_line_id = cur_line_id - 1

            if tau_trace_combined > tau_event and not disable_line_scattering:
                r_packet.next_line_id = cur_line_id
                r_packet.prev_line_id = cur_line_id - 1
                return (
                    distance_trace,
                    InteractionType.LINE,
                    delta_shell,
                )

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
            if distance_electron > 0.0:
                update_line_ids(
                    r_packet,
                    opacity_state.line_list_nu,
                    comov_nu_event,
                    line_id_step,
                )
            return (
                distance_electron,
                InteractionType.ESCATTERING,
                delta_shell,
            )
        interval_start = interval_end

    if distance_boundary > 0.0:
        update_line_ids(
            r_packet,
            opacity_state.line_list_nu,
            comov_nu_end,
            line_id_step,
        )
    return (
        distance_boundary,
        InteractionType.BOUNDARY,
        delta_shell,
    )
