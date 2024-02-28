from math import exp

from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)
from tardis.montecarlo.montecarlo_numba.numba_config import KB, H
from tardis.transport.frame_transformations import (
    calc_packet_energy,
    calc_packet_energy_full_relativity,
)


@njit(**njit_dict_no_parallel)
def update_base_estimators(
    r_packet, distance, estimator_state, comov_nu, comov_energy
):
    """
    Updating the estimators
    """
    estimator_state.j_estimator[r_packet.current_shell_id] += (
        comov_energy * distance
    )
    estimator_state.nu_bar_estimator[r_packet.current_shell_id] += (
        comov_energy * distance * comov_nu
    )


@njit(**njit_dict_no_parallel)
def update_bound_free_estimators(
    comov_nu,
    comov_energy,
    shell_id,
    distance,
    estimator_state,
    t_electron,
    x_sect_bfs,
    current_continua,
    bf_threshold_list_nu,
):
    """
    Update the estimators for bound-free processes.

    Parameters
    ----------
    comov_nu : float
    comov_energy : float
    shell_id : int
    distance : float
    numba_estimator : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    t_electron : float
        Electron temperature in the current cell.
    x_sect_bfs : numpy.ndarray, dtype float
        Photoionization cross-sections of all bound-free continua for
        which absorption is possible for frequency `comov_nu`.
    current_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `comov_nu`.
    bf_threshold_list_nu : numpy.ndarray, dtype float
        Threshold frequencies for photoionization sorted by decreasing frequency.
    """
    # TODO: Add full relativity mode
    boltzmann_factor = exp(-(H * comov_nu) / (KB * t_electron))
    for i, current_continuum in enumerate(current_continua):
        photo_ion_rate_estimator_increment = (
            comov_energy * distance * x_sect_bfs[i] / comov_nu
        )
        estimator_state.photo_ion_estimator[
            current_continuum, shell_id
        ] += photo_ion_rate_estimator_increment
        estimator_state.stim_recomb_estimator[current_continuum, shell_id] += (
            photo_ion_rate_estimator_increment * boltzmann_factor
        )
        estimator_state.photo_ion_estimator_statistics[
            current_continuum, shell_id
        ] += 1

        nu_th = bf_threshold_list_nu[current_continuum]
        bf_heating_estimator_increment = (
            comov_energy * distance * x_sect_bfs[i] * (1 - nu_th / comov_nu)
        )
        estimator_state.bf_heating_estimator[
            current_continuum, shell_id
        ] += bf_heating_estimator_increment
        estimator_state.stim_recomb_cooling_estimator[
            current_continuum, shell_id
        ] += (bf_heating_estimator_increment * boltzmann_factor)


@njit(**njit_dict_no_parallel)
def update_line_estimators(
    radfield_mc_estimators,
    r_packet,
    cur_line_id,
    distance_trace,
    time_explosion,
    enable_full_relativity,
):
    """
    Function to update the line estimators

    Parameters
    ----------
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    cur_line_id : int
    distance_trace : float
    time_explosion : float
    """
    if not enable_full_relativity:
        energy = calc_packet_energy(r_packet, distance_trace, time_explosion)
    else:
        energy = calc_packet_energy_full_relativity(r_packet)

    radfield_mc_estimators.j_blue_estimator[
        cur_line_id, r_packet.current_shell_id
    ] += (energy / r_packet.nu)
    radfield_mc_estimators.Edotlu_estimator[
        cur_line_id, r_packet.current_shell_id
    ] += energy
