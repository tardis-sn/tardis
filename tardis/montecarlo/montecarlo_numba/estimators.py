from math import exp

import numpy as np

from numba import njit, float64, int64
from numba.experimental import jitclass

from tardis.montecarlo import montecarlo_configuration as nc
from tardis.montecarlo.montecarlo_numba.numba_config import H, KB

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.transport.frame_transformations import (
    calc_packet_energy,
    calc_packet_energy_full_relativity,
)


def initialize_estimators(tau_sobolev_shape, gamma_shape):
    """
    Initializes the estimators used in the Monte Carlo simulation.

    Parameters
    ----------
    tau_sobolev_shape : tuple
        Shape of the array with the Sobolev optical depth.
    gamma_shape : tuple
        Shape of the array with the photoionization rate coefficients.

    Returns
    -------
    Estimators
        The initialized estimators.

    Examples
    --------
    >>> tau_sobolev_shape = (10, 20)
    >>> gamma_shape = (5, 5)
    >>> initialize_estimators(tau_sobolev_shape, gamma_shape)
    <Estimators object at 0x...>
    """

    j_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    nu_bar_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    j_blue_estimator = np.zeros(tau_sobolev_shape)
    Edotlu_estimator = np.zeros(tau_sobolev_shape)

    photo_ion_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_cooling_estimator = np.zeros(gamma_shape, dtype=np.float64)
    bf_heating_estimator = np.zeros(gamma_shape, dtype=np.float64)

    stim_recomb_cooling_estimator = np.zeros(gamma_shape, dtype=np.float64)

    photo_ion_estimator_statistics = np.zeros(gamma_shape, dtype=np.int64)
    return Estimators(
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        photo_ion_estimator_statistics,
    )


base_estimators_spec = [
    ("j_estimator", float64[:]),
    ("nu_bar_estimator", float64[:]),
    ("j_blue_estimator", float64[:, :]),
    ("Edotlu_estimator", float64[:, :]),
]

continuum_estimators_spec = [
    ("photo_ion_estimator", float64[:, :]),
    ("stim_recomb_estimator", float64[:, :]),
    ("bf_heating_estimator", float64[:, :]),
    ("stim_recomb_cooling_estimator", float64[:, :]),
    ("photo_ion_estimator_statistics", int64[:, :]),
]


@jitclass(base_estimators_spec + continuum_estimators_spec)
class Estimators(object):
    def __init__(
        self,
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        photo_ion_estimator_statistics,
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator
        self.photo_ion_estimator = photo_ion_estimator
        self.stim_recomb_estimator = stim_recomb_estimator
        self.bf_heating_estimator = bf_heating_estimator
        self.stim_recomb_cooling_estimator = stim_recomb_cooling_estimator
        self.photo_ion_estimator_statistics = photo_ion_estimator_statistics

    def increment(self, other):
        self.j_estimator += other.j_estimator
        self.nu_bar_estimator += other.nu_bar_estimator
        self.j_blue_estimator += other.j_blue_estimator
        self.Edotlu_estimator += other.Edotlu_estimator
        self.photo_ion_estimator += other.photo_ion_estimator
        self.stim_recomb_estimator += other.stim_recomb_estimator
        self.bf_heating_estimator += other.bf_heating_estimator
        self.stim_recomb_cooling_estimator += (
            other.stim_recomb_cooling_estimator
        )
        self.photo_ion_estimator_statistics += (
            other.photo_ion_estimator_statistics
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
    estimators, r_packet, cur_line_id, distance_trace, time_explosion
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

    """ Actual calculation - simplified below
    r_interaction = math.sqrt(r_packet.r**2 + distance_trace**2 +
                            2 * r_packet.r * distance_trace * r_packet.mu)
    mu_interaction = (r_packet.mu * r_packet.r + distance_trace) / r_interaction
    doppler_factor = 1.0 - mu_interaction * r_interaction /
    ( time_explosion * C)
    """

    if not nc.ENABLE_FULL_RELATIVITY:
        energy = calc_packet_energy(r_packet, distance_trace, time_explosion)
    else:
        energy = calc_packet_energy_full_relativity(r_packet)

    estimators.j_blue_estimator[cur_line_id, r_packet.current_shell_id] += (
        energy / r_packet.nu
    )
    estimators.Edotlu_estimator[
        cur_line_id, r_packet.current_shell_id
    ] += energy
