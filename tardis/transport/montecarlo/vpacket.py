import math

import numpy as np
from numba import float64, int64, njit
from numba.experimental import jitclass

from tardis.opacities.opacities import (
    chi_continuum_calculator,
)
from tardis.transport.frame_transformations import (
    angle_aberration_CMF_to_LF,
    angle_aberration_LF_to_CMF,
    get_doppler_factor,
)
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_line,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.numba_config import (
    C_SPEED_OF_LIGHT,
    SIGMA_THOMSON,
)
from tardis.transport.montecarlo.r_packet import (
    PacketStatus,
)
from tardis.transport.montecarlo.r_packet_transport import (
    move_packet_across_shell_boundary,
)

vpacket_spec = [
    ("r", float64),
    ("mu", float64),
    ("nu", float64),
    ("energy", float64),
    ("next_line_id", int64),
    ("current_shell_id", int64),
    ("status", int64),
    ("index", int64),
]


@jitclass(vpacket_spec)
class VPacket:
    def __init__(
        self,
        r,
        mu,
        nu,
        energy,
        current_shell_id,
        next_line_id,
        index=0,
    ):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = current_shell_id
        self.next_line_id = next_line_id
        self.status = PacketStatus.IN_PROCESS
        self.index = index


@njit(**njit_dict_no_parallel)
def trace_vpacket_within_shell(
    v_packet,
    numba_radial_1d_geometry,
    time_explosion,
    opacity_state,
    enable_full_relativity,
    continuum_processes_enabled,
):
    """
    Trace VPacket within one shell (relatively simple operation)
    """
    r_inner = numba_radial_1d_geometry.r_inner[v_packet.current_shell_id]
    r_outer = numba_radial_1d_geometry.r_outer[v_packet.current_shell_id]

    distance_boundary, delta_shell = calculate_distance_boundary(
        v_packet.r, v_packet.mu, r_inner, r_outer
    )
    # defining start for line interaction
    start_line_id = v_packet.next_line_id

    # e scattering initialization

    cur_electron_density = opacity_state.electron_density[
        v_packet.current_shell_id
    ]
    chi_e = cur_electron_density * SIGMA_THOMSON

    # Calculating doppler factor
    doppler_factor = get_doppler_factor(
        v_packet.r,
        v_packet.mu,
        time_explosion,
        enable_full_relativity,
    )

    comov_nu = v_packet.nu * doppler_factor

    if continuum_processes_enabled:
        (
            chi_bf_tot,
            chi_bf_contributions,
            current_continua,
            x_sect_bfs,
            chi_ff,
        ) = chi_continuum_calculator(
            opacity_state, comov_nu, v_packet.current_shell_id
        )
        chi_continuum = chi_e + chi_bf_tot + chi_ff

    else:
        chi_continuum = chi_e

    if enable_full_relativity:
        chi_continuum *= doppler_factor

    tau_continuum = chi_continuum * distance_boundary
    tau_trace_combined = tau_continuum

    cur_line_id = start_line_id

    for cur_line_id in range(start_line_id, len(opacity_state.line_list_nu)):
        # if tau_trace_combined > 10: ### FIXME ?????
        #    break

        nu_line = opacity_state.line_list_nu[cur_line_id]
        # TODO: Check if this is what the C code does

        tau_trace_line = opacity_state.tau_sobolev[
            cur_line_id, v_packet.current_shell_id
        ]

        is_last_line = cur_line_id == len(opacity_state.line_list_nu) - 1

        distance_trace_line = calculate_distance_line(
            v_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
            enable_full_relativity,
        )

        if distance_boundary <= distance_trace_line:
            break

        tau_trace_combined += tau_trace_line

    else:
        if cur_line_id == (len(opacity_state.line_list_nu) - 1):
            cur_line_id += 1
    v_packet.next_line_id = cur_line_id

    return tau_trace_combined, distance_boundary, delta_shell


@njit(**njit_dict_no_parallel)
def trace_vpacket(
    v_packet,
    numba_radial_1d_geometry,
    time_explosion,
    opacity_state,
    tau_russian,
    survival_probability,
    enable_full_relativity,
    continuum_processes_enabled,
):
    """
    Trace single vpacket.

    Parameters
    ----------
    v_packet
    time_explosion
    opacity_state

    Returns
    -------

    """
    tau_trace_combined = 0.0
    while True:
        (
            tau_trace_combined_shell,
            distance_boundary,
            delta_shell,
        ) = trace_vpacket_within_shell(
            v_packet,
            numba_radial_1d_geometry,
            time_explosion,
            opacity_state,
            enable_full_relativity,
            continuum_processes_enabled,
        )
        tau_trace_combined += tau_trace_combined_shell

        move_packet_across_shell_boundary(
            v_packet, delta_shell, len(numba_radial_1d_geometry.r_inner)
        )

        if tau_trace_combined > tau_russian:
            event_random = np.random.random()
            if event_random > survival_probability:
                v_packet.energy = 0.0
                v_packet.status = PacketStatus.EMITTED
            else:
                v_packet.energy = (
                    v_packet.energy
                    / survival_probability
                    * math.exp(-tau_trace_combined)
                )
                tau_trace_combined = 0.0

        # Moving the v_packet
        new_r = math.sqrt(
            v_packet.r * v_packet.r
            + distance_boundary * distance_boundary
            + 2.0 * v_packet.r * distance_boundary * v_packet.mu
        )
        v_packet.mu = (v_packet.mu * v_packet.r + distance_boundary) / new_r
        v_packet.r = new_r

        if v_packet.status == PacketStatus.EMITTED:
            break
    return tau_trace_combined


@njit(**njit_dict_no_parallel)
def trace_vpacket_volley(
    r_packet,
    vpacket_collection,
    numba_radial_1d_geometry,
    time_explosion,
    opacity_state,
    enable_full_relativity,
    tau_russian,
    survival_probability,
    continuum_processes_enabled,
):
    """
    Shoot a volley of vpackets (the vpacket collection specifies how many)
    from the current position of the rpacket.

    Parameters
    ----------
    r_packet : [type]
        [description]
    vpacket_collection : [type]
        [description]
    numba_radial_1d_geometry : [type]
        [description]
    time_explosion : [type]
        [description]
    opacity_state : [type]
        [description]
    """
    if (r_packet.nu < vpacket_collection.v_packet_spawn_start_frequency) or (
        r_packet.nu > vpacket_collection.v_packet_spawn_end_frequency
    ):
        return

    no_of_vpackets = vpacket_collection.number_of_vpackets
    if no_of_vpackets == 0:
        return

    ### TODO theoretical check for r_packet nu within vpackets bins - is done somewhere else I think
    if (
        r_packet.r > numba_radial_1d_geometry.r_inner[0]
    ):  # not on inner_boundary
        r_inner_over_r = numba_radial_1d_geometry.r_inner[0] / r_packet.r
        mu_min = -math.sqrt(1 - r_inner_over_r * r_inner_over_r)
        v_packet_on_inner_boundary = False
        if enable_full_relativity:
            mu_min = angle_aberration_LF_to_CMF(
                r_packet, time_explosion, mu_min
            )
    else:
        v_packet_on_inner_boundary = True
        mu_min = 0.0

        if enable_full_relativity:
            inv_c = 1 / C_SPEED_OF_LIGHT
            inv_t = 1 / time_explosion
            beta_inner = numba_radial_1d_geometry.r_inner[0] * inv_t * inv_c

    mu_bin = (1.0 - mu_min) / no_of_vpackets
    r_packet_doppler_factor = get_doppler_factor(
        r_packet.r,
        r_packet.mu,
        time_explosion,
        enable_full_relativity,
    )
    for i in range(no_of_vpackets):
        v_packet_mu = mu_min + i * mu_bin + np.random.random() * mu_bin

        if v_packet_on_inner_boundary:  # The weights are described in K&S 2014
            if not enable_full_relativity:
                weight = 2 * v_packet_mu / no_of_vpackets
            else:
                weight = (
                    2
                    * (v_packet_mu + beta_inner)
                    / (2 * beta_inner + 1)
                    / no_of_vpackets
                )

        else:
            weight = (1 - mu_min) / (2 * no_of_vpackets)

        # C code: next line, angle_aberration_CMF_to_LF( & virt_packet, storage);
        if enable_full_relativity:
            v_packet_mu = angle_aberration_CMF_to_LF(
                r_packet, time_explosion, v_packet_mu
            )
        v_packet_doppler_factor = get_doppler_factor(
            r_packet.r,
            v_packet_mu,
            time_explosion,
            enable_full_relativity,
        )

        # transform between r_packet mu and v_packet_mu

        doppler_factor_ratio = r_packet_doppler_factor / v_packet_doppler_factor

        v_packet_nu = r_packet.nu * doppler_factor_ratio
        v_packet_energy = r_packet.energy * weight * doppler_factor_ratio

        # TODO: Make sure we have a new continuum object for each vpacket
        # comov_nu = v_packet_nu * v_packet_doppler_factor
        # continuum.calculate(comov_nu, r_packet.current_shell_id)

        v_packet = VPacket(
            r_packet.r,
            v_packet_mu,
            v_packet_nu,
            v_packet_energy,
            r_packet.current_shell_id,
            r_packet.next_line_id,
            i,
        )

        tau_vpacket = trace_vpacket(
            v_packet,
            numba_radial_1d_geometry,
            time_explosion,
            opacity_state,
            tau_russian,
            survival_probability,
            enable_full_relativity,
            continuum_processes_enabled,
        )

        v_packet.energy *= math.exp(-tau_vpacket)

        vpacket_collection.add_packet(
            v_packet.nu,
            v_packet.energy,
            v_packet_mu,
            r_packet.r,
            r_packet.last_interaction_in_nu,
            r_packet.last_interaction_type,
            r_packet.last_line_interaction_in_id,
            r_packet.last_line_interaction_out_id,
            r_packet.last_line_interaction_shell_id,
        )
