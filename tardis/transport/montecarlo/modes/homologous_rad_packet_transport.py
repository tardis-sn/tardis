"""Shared homologous-mode radiative packet tracing."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from tardis.transport.frame_transformations import get_doppler_factor
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_line,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_line,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
)

if TYPE_CHECKING:
    from tardis.model.geometry.radial1d_nonhomologous import (
        NumbaNonhomologousRadial1DGeometry,
    )
    from tardis.transport.montecarlo.estimators.estimators_line import (
        EstimatorsLine,
    )
    from tardis.transport.montecarlo.packets.radiative_packet import RPacket


@njit(**njit_dict_no_parallel)
def trace_packet(
    r_packet: RPacket,
    numba_radial_1d_geometry: NumbaNonhomologousRadial1DGeometry,
    time_explosion: float,
    opacity_state,
    estimators_line: EstimatorsLine,
    continuous_opacity: float,
    escat_prob: float,
    continuum_process_enabled: bool,
    enable_full_relativity: bool,
    disable_line_scattering: bool,
) -> tuple[float, InteractionType, int]:
    """
    Trace a homologous radiative packet until an interaction occurs.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet being transported.
    numba_radial_1d_geometry : NumbaNonhomologousRadial1DGeometry
        Radial 1D geometry of the model.
    time_explosion : float
        Time since explosion in seconds.
    opacity_state
        Opacity state containing line frequencies and Sobolev optical depths.
    estimators_line : EstimatorsLine
        Line-level radiation field estimators.
    continuous_opacity : float
        Continuous interaction opacity. In classic mode this is electron
        scattering opacity; in IIP mode this is total continuum opacity.
    escat_prob : float
        Probability that a continuous interaction is electron scattering.
        Ignored when continuum processes are disabled.
    continuum_process_enabled : bool
        Whether continuous interactions may produce continuum-process events.
    enable_full_relativity : bool
        Whether full relativistic calculations are enabled.
    disable_line_scattering : bool
        Whether line scattering is disabled.

    Returns
    -------
    tuple[float, InteractionType, int]
        Distance, interaction type, and shell-index change.
    """
    r_inner = numba_radial_1d_geometry.r_inner[r_packet.current_shell_id]
    r_outer = numba_radial_1d_geometry.r_outer[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    start_line_id = r_packet.next_line_id
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    velocity = numba_radial_1d_geometry.get_velocity(
        r_packet.r, r_packet.current_shell_id
    )
    doppler_factor = get_doppler_factor(
        velocity,
        r_packet.mu,
        enable_full_relativity,
    )
    comov_nu = r_packet.nu * doppler_factor

    distance_continuous = tau_event / continuous_opacity
    cur_line_id = start_line_id
    last_line_id = len(opacity_state.line_list_nu) - 1
    for cur_line_id in range(start_line_id, len(opacity_state.line_list_nu)):
        nu_line = opacity_state.line_list_nu[cur_line_id]
        tau_trace_line = opacity_state.tau_sobolev[
            cur_line_id, r_packet.current_shell_id
        ]
        tau_trace_line_combined += tau_trace_line

        is_last_line = cur_line_id == last_line_id
        distance_trace = calculate_distance_line(
            r_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
            enable_full_relativity,
        )

        tau_trace_continuous = continuous_opacity * distance_trace
        tau_trace_combined = (
            tau_trace_line_combined + tau_trace_continuous
        )
        distance = min(
            distance_trace, distance_boundary, distance_continuous
        )

        if distance_trace != 0:
            if distance == distance_boundary:
                interaction_type = InteractionType.BOUNDARY
                r_packet.next_line_id = cur_line_id
                break
            if distance == distance_continuous:
                if continuum_process_enabled:
                    zrand = np.random.random()
                    if zrand < escat_prob:
                        interaction_type = InteractionType.ESCATTERING
                    else:
                        interaction_type = InteractionType.CONTINUUM_PROCESS
                else:
                    interaction_type = InteractionType.ESCATTERING
                r_packet.next_line_id = cur_line_id
                break

        update_estimators_line(
            estimators_line,
            r_packet,
            cur_line_id,
            distance_trace,
            time_explosion,
            enable_full_relativity,
        )

        if tau_trace_combined > tau_event and not disable_line_scattering:
            interaction_type = InteractionType.LINE
            r_packet.next_line_id = cur_line_id
            distance = distance_trace
            break

        distance_continuous = (
            tau_event - tau_trace_line_combined
        ) / continuous_opacity

    else:
        if cur_line_id == last_line_id:
            cur_line_id += 1
        if distance_continuous < distance_boundary:
            distance = distance_continuous
            if continuum_process_enabled:
                zrand = np.random.random()
                if zrand < escat_prob:
                    interaction_type = InteractionType.ESCATTERING
                else:
                    interaction_type = InteractionType.CONTINUUM_PROCESS
            else:
                interaction_type = InteractionType.ESCATTERING
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    return distance, interaction_type, delta_shell
