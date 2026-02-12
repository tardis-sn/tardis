from math import exp

from numba import njit

from tardis.transport.frame_transformations import (
    calc_packet_energy,
    calc_packet_energy_full_relativity,
)
from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)
from tardis.transport.montecarlo.configuration.constants import KB, H
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)


@njit(**njit_dict_no_parallel)
def update_estimators_bulk(
    r_packet,
    distance: float,
    estimators_bulk: EstimatorsBulk,
    comov_nu: float,
    comov_energy: float,
) -> None:
    """
    Update cell-level bulk radiation field estimators.

    Parameters
    ----------
    r_packet
        The radiative packet being transported.
    distance
        Distance traveled by the packet.
    estimators_bulk
        Cell-level bulk radiation field estimators.
    comov_nu
        Comoving frame frequency.
    comov_energy
        Comoving frame energy.
    """
    estimators_bulk.mean_intensity_total[r_packet.current_shell_id] += (
        comov_energy * distance
    )
    estimators_bulk.mean_frequency[r_packet.current_shell_id] += (
        comov_energy * distance * comov_nu
    )


@njit(**njit_dict_no_parallel)
def update_estimators_bound_free(
    comov_nu: float,
    comov_energy: float,
    shell_id: int,
    distance: float,
    estimators_continuum: EstimatorsContinuum,
    t_electron: float,
    x_sect_bfs,
    current_continua,
    bf_threshold_list_nu,
    chi_ff: float,
) -> None:
    """
    Update the estimators for bound-free processes in place by thread.

    Parameters
    ----------
    comov_nu
        Comoving frame frequency.
    comov_energy
        Comoving frame energy.
    shell_id
        Current cell index.
    distance
        Distance traveled by the packet in the current shell.
    estimators_continuum
        Continuum interaction estimators.
    t_electron
        Electron temperature in the current cell.
    x_sect_bfs
        Photoionization cross-sections of all bound-free continua for
        which absorption is possible for frequency `comov_nu`.
    current_continua
        Continuum ids for which absorption is possible for frequency `comov_nu`.
    bf_threshold_list_nu
        Threshold frequencies for photoionization sorted by decreasing frequency.
    chi_ff
        Free-free opacity coefficient in the current cell.
    """
    # TODO: Add full relativity mode
    boltzmann_factor = exp(-(H * comov_nu) / (KB * t_electron))
    estimators_continuum.ff_heating_estimator[shell_id] += (
        comov_energy * distance * chi_ff
    )
    for i, current_continuum in enumerate(current_continua):
        photo_ion_rate_estimator_increment = (
            comov_energy * distance * x_sect_bfs[i] / comov_nu
        )
        estimators_continuum.photo_ion_estimator[
            current_continuum, shell_id
        ] += photo_ion_rate_estimator_increment
        estimators_continuum.stim_recomb_estimator[
            current_continuum, shell_id
        ] += photo_ion_rate_estimator_increment * boltzmann_factor
        estimators_continuum.photo_ion_estimator_statistics[
            current_continuum, shell_id
        ] += 1

        nu_th = bf_threshold_list_nu[current_continuum]
        bf_heating_estimator_increment = (
            comov_energy * distance * x_sect_bfs[i] * (1 - nu_th / comov_nu)
        )
        estimators_continuum.bf_heating_estimator[
            current_continuum, shell_id
        ] += bf_heating_estimator_increment
        estimators_continuum.stim_recomb_cooling_estimator[
            current_continuum, shell_id
        ] += bf_heating_estimator_increment * boltzmann_factor


@njit(**njit_dict_no_parallel)
def update_estimators_line(
    estimators_line: EstimatorsLine,
    r_packet,
    cur_line_id: int,
    distance_trace: float,
    time_explosion: float,
    enable_full_relativity: bool,
) -> None:
    """
    Update line-level radiation field estimators in place by thread.

    Parameters
    ----------
    estimators_line
        Line-level radiation field estimators.
    r_packet
        The radiative packet being transported.
    cur_line_id
        Current line index.
    distance_trace
        Distance traced by the packet.
    time_explosion
        Time since explosion in seconds.
    enable_full_relativity
        Flag to enable full relativistic treatment.
    """
    if not enable_full_relativity:
        energy = calc_packet_energy(r_packet, distance_trace, time_explosion)
    else:
        energy = calc_packet_energy_full_relativity(r_packet)

    estimators_line.mean_intensity_blueward[
        cur_line_id, r_packet.current_shell_id
    ] += energy / r_packet.nu
    estimators_line.energy_deposition_line_rate[
        cur_line_id, r_packet.current_shell_id
    ] += energy
