from numba import njit
from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType)

from tardis.montecarlo import montecarlo_configuration as montecarlo_configuration
from tardis.montecarlo.montecarlo_numba.r_packet import (
    get_doppler_factor, get_inverse_doppler_factor, get_random_mu,
    angle_aberration_CMF_to_LF, test_for_close_line)
from tardis.montecarlo.montecarlo_numba.macro_atom import macro_atom

@njit(**njit_dict)
def thomson_scatter(r_packet, time_explosion):
    """
    Thomson scattering — no longer line scattering
    2) get the doppler factor at that position with the old angle
    3) convert the current energy and nu into the comoving
        frame with the old mu
    4) Scatter and draw new mu - update mu
    5) Transform the comoving energy and nu back using the new mu

    Parameters
    ----------
    distance : [type]
        [description]
    """
    old_doppler_factor = get_doppler_factor(
        r_packet.r,
        r_packet.mu,
        time_explosion)
    comov_nu = r_packet.nu * old_doppler_factor
    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion)

    r_packet.nu = comov_nu * inverse_new_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet,
            time_explosion,
            r_packet.mu
        )


@njit(**njit_dict)
def line_scatter(r_packet, time_explosion, line_interaction_type, numba_plasma):
    #increment_j_blue_estimator(packet, storage, distance, line2d_idx);
    #increment_Edotlu_estimator(packet, storage, distance, line2d_idx);

    # do_electron_scatter = False
    # general_scatter(r_packet, time_explosion, do_electron_scatter)
    # update last_interaction

    old_doppler_factor = get_doppler_factor(r_packet.r,
                                        r_packet.mu,
                                        time_explosion)
    r_packet.mu = get_random_mu()

    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion)

    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor

    if line_interaction_type == LineInteractionType.SCATTER:
        line_emission(r_packet, r_packet.next_line_id,
                      time_explosion, numba_plasma)
    else: # includes both macro atom and downbranch - encoded in the transition probabilities
        emission_line_id = macro_atom(r_packet, numba_plasma)
        line_emission(r_packet, emission_line_id, time_explosion,
                      numba_plasma)

@njit(**njit_dict)
def line_emission(r_packet, emission_line_id, time_explosion,
                  numba_plasma):
    if emission_line_id != r_packet.next_line_id:
        pass
    inverse_doppler_factor = get_inverse_doppler_factor(r_packet.r, r_packet.mu,
                                        time_explosion)
    r_packet.nu = numba_plasma.line_list_nu[
                      emission_line_id] * inverse_doppler_factor
    r_packet.next_line_id = emission_line_id + 1
    nu_line = numba_plasma.line_list_nu[emission_line_id]

    if emission_line_id != (len(numba_plasma.line_list_nu) - 1):
        test_for_close_line(r_packet, emission_line_id + 1, nu_line, numba_plasma)

    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet,
            time_explosion,
            r_packet.mu
            )




