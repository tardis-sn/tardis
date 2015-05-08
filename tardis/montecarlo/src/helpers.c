#include "cmontecarlo.h"

rk_state mt_state;


INLINE double
rpacket_doppler_factor (rpacket_t * packet, storage_model_t * storage)
{
  return 1.0 -
    rpacket_get_mu (packet) * rpacket_get_r (packet) *
    storage->inverse_time_explosion * INVERSE_C;
}


void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
			 double distance)
{
  double comov_energy = 0.0;
  int64_t emission_line_id = 0;
  double old_doppler_factor = 0.0;
  double inverse_doppler_factor = 0.0;
  double tau_line = 0.0;
  double tau_electron = 0.0;
  double tau_combined = 0.0;
  bool virtual_close_line = false;
  int64_t j_blue_idx = -1;
  if (rpacket_get_virtual_packet (packet) == 0)
    {
      j_blue_idx =
	rpacket_get_current_shell_id (packet) *
	storage->line_lists_j_blues_nd + rpacket_get_next_line_id (packet);
      increment_j_blue_estimator (packet, storage, distance, j_blue_idx);
    }
  tau_line =
    storage->line_lists_tau_sobolevs[rpacket_get_current_shell_id (packet) *
				     storage->line_lists_tau_sobolevs_nd +
				     rpacket_get_next_line_id (packet)];
  tau_electron =
    storage->sigma_thomson *
    storage->electron_densities[rpacket_get_current_shell_id (packet)] *
    distance;
  tau_combined = tau_line + tau_electron;
  rpacket_set_next_line_id (packet, rpacket_get_next_line_id (packet) + 1);
  if (rpacket_get_next_line_id (packet) == storage->no_of_lines)
    {
      rpacket_set_last_line (packet, true);
    }
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) + tau_line);
    }
  else if (rpacket_get_tau_event (packet) < tau_combined)
    {
      old_doppler_factor = move_packet (packet, storage, distance);
      rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
      inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      comov_energy = rpacket_get_energy (packet) * old_doppler_factor;
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      storage->last_line_interaction_in_id[storage->current_packet_id] =
	rpacket_get_next_line_id (packet) - 1;
      storage->last_line_interaction_shell_id[storage->current_packet_id] =
	rpacket_get_current_shell_id (packet);
      storage->last_interaction_type[storage->current_packet_id] = 2;
      if (storage->line_interaction_id == 0)
	{
	  emission_line_id = rpacket_get_next_line_id (packet) - 1;
	}
      else if (storage->line_interaction_id >= 1)
	{
	  emission_line_id = macro_atom (packet, storage);
	}
      storage->last_line_interaction_out_id[storage->current_packet_id] =
	emission_line_id;
      rpacket_set_nu (packet,
		      storage->line_list_nu[emission_line_id] *
		      inverse_doppler_factor);
      rpacket_set_nu_line (packet, storage->line_list_nu[emission_line_id]);
      rpacket_set_next_line_id (packet, emission_line_id + 1);
      rpacket_reset_tau_event (packet);
      rpacket_set_recently_crossed_boundary (packet, 0);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  virtual_close_line = false;
	  if (!rpacket_get_last_line (packet) &&
	      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
		    rpacket_get_nu_line (packet)) /
	      rpacket_get_nu_line (packet) < 1e-7)
	    {
	      virtual_close_line = true;
	    }
	  // QUESTIONABLE!!!
	  bool old_close_line = rpacket_get_close_line (packet);
	  rpacket_set_close_line (packet, virtual_close_line);
	  montecarlo_one_packet (storage, packet, 1);
	  rpacket_set_close_line (packet, old_close_line);
	  virtual_close_line = false;
	}
    }
  else
    {
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) - tau_line);
    }
  if (!rpacket_get_last_line (packet) &&
      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
	    rpacket_get_nu_line (packet)) / rpacket_get_nu_line (packet) <
      1e-7)
    {
      rpacket_set_close_line (packet, true);
    }
}

void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
			    double distance)
{
  double comov_energy, doppler_factor, comov_nu, inverse_doppler_factor;
  doppler_factor = move_packet (packet, storage, distance);
  comov_nu = rpacket_get_nu (packet) * doppler_factor;
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  rpacket_set_mu (packet, 2.0 * rk_double (&mt_state) - 1.0);
  inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
  rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
  rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
  rpacket_reset_tau_event (packet);
  rpacket_set_recently_crossed_boundary (packet, 0);
  storage->last_interaction_type[storage->current_packet_id] = 1;
  if (rpacket_get_virtual_packet_flag (packet) > 0)
    {
      montecarlo_one_packet (storage, packet, 1);
    }
}

void
move_packet_across_shell_boundary (rpacket_t * packet,
				   storage_model_t * storage, double distance)
{
  double comov_energy, doppler_factor, comov_nu, inverse_doppler_factor;
  move_packet (packet, storage, distance);
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      double delta_tau_event = distance *
	storage->electron_densities[rpacket_get_current_shell_id (packet)] *
	storage->sigma_thomson;
      rpacket_set_tau_event (packet,
			     rpacket_get_tau_event (packet) +
			     delta_tau_event);
    }
  else
    {
      rpacket_reset_tau_event (packet);
    }
  if ((rpacket_get_current_shell_id (packet) < storage->no_of_shells - 1
       && rpacket_get_next_shell_id (packet) == 1)
      || (rpacket_get_current_shell_id (packet) > 0
	  && rpacket_get_next_shell_id (packet) == -1))
    {
      rpacket_set_current_shell_id (packet,
				    rpacket_get_current_shell_id (packet) +
				    rpacket_get_next_shell_id (packet));
      rpacket_set_recently_crossed_boundary (packet,
					     rpacket_get_next_shell_id
					     (packet));
    }
  else if (rpacket_get_next_shell_id (packet) == 1)
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
    }
  else if ((storage->reflective_inner_boundary == 0) ||
	   (rk_double (&mt_state) > storage->inner_boundary_albedo))
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
    }
  else
    {
      doppler_factor = rpacket_doppler_factor (packet, storage);
      comov_nu = rpacket_get_nu (packet) * doppler_factor;
      comov_energy = rpacket_get_energy (packet) * doppler_factor;
      rpacket_set_mu (packet, rk_double (&mt_state));
      inverse_doppler_factor = 1.0 / rpacket_doppler_factor (packet, storage);
      rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      rpacket_set_recently_crossed_boundary (packet, 1);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
	{
	  montecarlo_one_packet (storage, packet, -2);
	}
    }
}

INLINE montecarlo_event_handler_t
get_event_handler (rpacket_t * packet, storage_model_t * storage,
		   double *distance)
{
  double d_boundary, d_electron, d_line;
  montecarlo_compute_distances (packet, storage);
  d_boundary = rpacket_get_d_boundary (packet);
  d_electron = rpacket_get_d_electron (packet);
  d_line = rpacket_get_d_line (packet);
  montecarlo_event_handler_t handler;
  if (d_line <= d_boundary && d_line <= d_electron)
    {
      *distance = d_line;
      handler = &montecarlo_line_scatter;
    }
  else if (d_boundary <= d_electron)
    {
      *distance = d_boundary;
      handler = &move_packet_across_shell_boundary;
    }
  else
    {
      *distance = d_electron;
      handler = &montecarlo_thomson_scatter;
    }
  return handler;
}

