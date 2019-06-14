from enum import IntEnum

from numba import njit

from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis.montecarlo.montecarlo_numba.rpacket import get_doppler_factor, \
    get_random_mu

class LineInteractionType(IntEnum):
    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2

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
    old_doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    comov_energy = r_packet.energy * old_doppler_factor
    comov_nu = r_packet.nu * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = 1. / get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion)
    r_packet.energy = comov_energy * inverse_new_doppler_factor
    r_packet.nu = comov_nu * inverse_new_doppler_factor

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

def line_scatter(r_packet, time_explosion, line_interaction_type, numba_plasma):
    #increment_j_blue_estimator(packet, storage, distance, line2d_idx);
    #increment_Edotlu_estimator(packet, storage, distance, line2d_idx);

    general_scatter(r_packet, time_explosion)
    # update last_interaction

    if line_interaction_type == LineInteractionType.SCATTER:
        line_emission(r_packet, r_packet.next_line_id, time_explosion, numba_plasma)
    else:
        pass
"""        
        line_emission()

    if (storage->line_interaction_id == 0)
    {
        line_emission(packet, storage, next_line_id, mt_state);
    }
    else if (storage->line_interaction_id >= 1)
        {
        rpacket_set_macro_atom_activation_level(packet,
                                                    storage->
        rpacket_set_macro_atom_activation_level(packet,
                                                storage->line2macro_level_upper[
            next_line_id]);
        macro_atom(packet, storage, mt_state);
"""

def line_emission(r_packet, emission_line_id, time_explosion,
                  numba_plasma):
    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        time_explosion)
    r_packet.nu = numba_plasma.line_list_nu[
                      emission_line_id] / doppler_factor
    r_packet.next_line_id = emission_line_id + 1

"""
void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
                         double distance, rk_state *mt_state)
{
  uint64_t next_line_id = rpacket_get_next_line_id (packet);
  uint64_t line2d_idx = next_line_id +
    storage->no_of_lines * rpacket_get_current_shell_id (packet);
  if (rpacket_get_virtual_packet (packet) == 0)
    {
    }
  double tau_line =
    storage->line_lists_tau_sobolevs[line2d_idx];
  double tau_continuum = rpacket_get_chi_continuum(packet) * distance;
  double tau_combined = tau_line + tau_continuum;
  //rpacket_set_next_line_id (packet, rpacket_get_next_line_id (packet) + 1);

  if (next_line_id + 1 == storage->no_of_lines)
    {
      rpacket_set_last_line (packet, true);
    }
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      rpacket_set_tau_event (packet,
                             rpacket_get_tau_event (packet) + tau_line);
      rpacket_set_next_line_id (packet, next_line_id + 1);
      test_for_close_line (packet, storage);
    }
  else if (rpacket_get_tau_event (packet) < tau_combined)
    { // Line absorption occurs
      move_packet (packet, storage, distance);
      double old_doppler_factor = rpacket_doppler_factor (packet, storage);
      rpacket_set_mu (packet, 2.0 * rk_double (mt_state) - 1.0);
      double inverse_doppler_factor = rpacket_inverse_doppler_factor (packet, storage);
      double comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_interaction_in_nu[rpacket_get_id (packet)] =
        rpacket_get_nu (packet);
      storage->last_line_interaction_in_id[rpacket_get_id (packet)] =
        next_line_id;
      storage->last_line_interaction_shell_id[rpacket_get_id (packet)] =
        rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[rpacket_get_id (packet)] = 2;
      if (storage->line_interaction_id == 0)
        {
          line_emission (packet, storage, next_line_id, mt_state);
        }
      else if (storage->line_interaction_id >= 1)
        {
          rpacket_set_macro_atom_activation_level (packet,
                                                   storage->line2macro_level_upper[next_line_id]);
          macro_atom (packet, storage, mt_state);
        }
    }
  else
    { // Packet passes line without interacting
      rpacket_set_tau_event (packet,
                             rpacket_get_tau_event (packet) - tau_line);
      rpacket_set_next_line_id (packet, next_line_id + 1);
      packet->compute_chi_bf = false;
      test_for_close_line (packet, storage);
    }
}
"""

