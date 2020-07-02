from numba import njit
from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType)

from tardis.montecarlo import montecarlo_configuration as montecarlo_configuration
from tardis.montecarlo.montecarlo_numba.r_packet import (
    get_doppler_factor, get_inverse_doppler_factor, get_random_mu,
    angle_aberration_CMF_to_LF)
from tardis.montecarlo.montecarlo_numba.macro_atom import macro_atom


@njit(**njit_dict)
def general_scatter(r_packet, time_explosion):
    """
    Thomson as well as line scattering
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
    comov_energy = r_packet.energy * old_doppler_factor
    comov_nu = r_packet.nu * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion)
    r_packet.energy = comov_energy * inverse_new_doppler_factor
    r_packet.nu = comov_nu * inverse_new_doppler_factor
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet,
            time_explosion,
            r_packet.mu
        )


"""
void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
                            double distance, rk_state *mt_state)
{
  move_packet (packet, storage, distance);
  double doppler_factor = rpacket_doppler_factor (packet, storage);
  double comov_nu = rpacket_get_nu (packet) * doppler_factor;
  double comov_energy = rpacket_get_energy (packet) * doppler_factor;
  rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
  double inverse_doppler_factor = rpacket_inverse_doppler_factor (packet, storage);
  rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  rpacket_reset_tau_event (packet, mt_state);
  storage->last_interaction_type[rpacket_get_id (packet)] = 1;

  angle_aberration_CMF_to_LF (packet, storage);

  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      create_vpacket (storage, packet, mt_state);
    }
}

"""

@njit(**njit_dict)
def line_scatter(r_packet, time_explosion, line_interaction_type, numba_plasma):
    #increment_j_blue_estimator(packet, storage, distance, line2d_idx);
    #increment_Edotlu_estimator(packet, storage, distance, line2d_idx);

    general_scatter(r_packet, time_explosion)
    # update last_interaction

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
    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        time_explosion)
    r_packet.nu = numba_plasma.line_list_nu[
                      emission_line_id] / doppler_factor
    r_packet.next_line_id = emission_line_id + 1
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet,
            time_explosion,
            r_packet.mu
            )




