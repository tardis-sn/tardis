from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.montecarlo_numba.frame_transformations import (
    calc_packet_energy,
    calc_packet_energy_full_relativity,
)

from tardis.montecarlo.montecarlo_numba.numba_config import (
    ENABLE_FULL_RELATIVITY,
)


@njit(**njit_dict_no_parallel)
def set_estimators(r_packet, distance, numba_estimator, comov_nu, comov_energy):
    """
    Updating the estimators
    """
    numba_estimator.j_estimator[r_packet.current_shell_id] += (
        comov_energy * distance
    )
    numba_estimator.nu_bar_estimator[r_packet.current_shell_id] += (
        comov_energy * distance * comov_nu
    )


@njit(**njit_dict_no_parallel)
def set_estimators_full_relativity(
    r_packet, distance, numba_estimator, comov_nu, comov_energy, doppler_factor
):
    numba_estimator.j_estimator[r_packet.current_shell_id] += (
        comov_energy * distance * doppler_factor
    )
    numba_estimator.nu_bar_estimator[r_packet.current_shell_id] += (
        comov_energy * distance * comov_nu * doppler_factor
    )


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

    if not ENABLE_FULL_RELATIVITY:
        energy = calc_packet_energy(r_packet, distance_trace, time_explosion)
    else:
        energy = calc_packet_energy_full_relativity(r_packet)

    estimators.j_blue_estimator[cur_line_id, r_packet.current_shell_id] += (
        energy / r_packet.nu
    )
    estimators.Edotlu_estimator[
        cur_line_id, r_packet.current_shell_id
    ] += energy
