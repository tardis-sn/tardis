#include "cmontecarlo.h"
#include ""

rk_state mt_state;

void
initialize_random_kit (unsigned long seed)
{
  rk_seed (seed, &mt_state);
}

INLINE int64_t
macro_atom (rpacket_t * packet, storage_model_t * storage)
{
  int emit = 0, i = 0;
  double p, event_random;
  int activate_level =
    storage->line2macro_level_upper[rpacket_get_next_line_id (packet) - 1];
  while (emit != -1)
    {
      event_random = rk_double (&mt_state);
      i = storage->macro_block_references[activate_level] - 1;
      p = 0.0;
      do
	{
	  p +=
	    storage->
	    transition_probabilities[rpacket_get_current_shell_id (packet) *
				     storage->transition_probabilities_nd +
				     (++i)];
	}
      while (p <= event_random);
      emit = storage->transition_type[i];
      activate_level = storage->destination_level_id[i];
    }
  return storage->transition_line_id[i];
}

INLINE double
move_packet (rpacket_t * packet, storage_model_t * storage, double distance)
{
  double new_r, doppler_factor, comov_energy, comov_nu;
  doppler_factor = rpacket_doppler_factor (packet, storage);
  if (distance > 0.0)
    {
      double r = rpacket_get_r (packet);
      new_r =
	sqrt (r * r + distance * distance +
	      2.0 * r * distance * rpacket_get_mu (packet));
      rpacket_set_mu (packet,
		      (rpacket_get_mu (packet) * r + distance) / new_r);
      rpacket_set_r (packet, new_r);
      if (rpacket_get_virtual_packet (packet) <= 0)
	{
	  comov_energy = rpacket_get_energy (packet) * doppler_factor;
	  comov_nu = rpacket_get_nu (packet) * doppler_factor;
	  storage->js[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance;
	  storage->nubars[rpacket_get_current_shell_id (packet)] +=
	    comov_energy * distance * comov_nu;
	}
    }
  return doppler_factor;
}

INLINE void
increment_j_blue_estimator (rpacket_t * packet, storage_model_t * storage,
			    double d_line, int64_t j_blue_idx)
{
  double comov_energy, r_interaction, mu_interaction, doppler_factor;
  double r = rpacket_get_r (packet);
  r_interaction =
    sqrt (r * r + d_line * d_line +
	  2.0 * r * d_line * rpacket_get_mu (packet));
  mu_interaction = (rpacket_get_mu (packet) * r + d_line) / r_interaction;
  doppler_factor = 1.0 - mu_interaction * r_interaction *
    storage->inverse_time_explosion * INVERSE_C;
  comov_energy = rpacket_get_energy (packet) * doppler_factor;
  storage->line_lists_j_blues[j_blue_idx] +=
    comov_energy / rpacket_get_nu (packet);
}

int64_t
montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
		       int64_t virtual_mode)
{
  int64_t i;
  rpacket_t virt_packet;
  double mu_bin;
  double mu_min;
  double doppler_factor_ratio;
  double weight;
  int64_t virt_id_nu;
  int64_t reabsorbed;
  if (virtual_mode == 0)
    {
      reabsorbed = montecarlo_one_packet_loop (storage, packet, 0);
    }
  else
    {
      for (i = 0; i < rpacket_get_virtual_packet_flag (packet); i++)
	{
	  memcpy ((void *) &virt_packet, (void *) packet, sizeof (rpacket_t));
	  if (virt_packet.r > storage->r_inner[0])
	    {
	      mu_min =
		-1.0 * sqrt (1.0 -
			     (storage->r_inner[0] / virt_packet.r) *
			     (storage->r_inner[0] / virt_packet.r));
	    }
	  else
	    {
	      mu_min = 0.0;
	    }
	  mu_bin = (1.0 - mu_min) / rpacket_get_virtual_packet_flag (packet);
	  virt_packet.mu = mu_min + (i + rk_double (&mt_state)) * mu_bin;
	  switch (virtual_mode)
	    {
	    case -2:
	      weight = 1.0 / rpacket_get_virtual_packet_flag (packet);
	      break;
	    case -1:
	      weight =
		2.0 * virt_packet.mu /
		rpacket_get_virtual_packet_flag (packet);
	      break;
	    case 1:
	      weight =
		(1.0 -
		 mu_min) / 2.0 / rpacket_get_virtual_packet_flag (packet);
	      break;
	    default:
	      fprintf (stderr, "Something has gone horribly wrong!\n");
	    }
	  doppler_factor_ratio =
	    rpacket_doppler_factor (packet, storage) /
	    rpacket_doppler_factor (&virt_packet, storage);
	  virt_packet.energy =
	    rpacket_get_energy (packet) * doppler_factor_ratio;
	  virt_packet.nu = rpacket_get_nu (packet) * doppler_factor_ratio;
	  reabsorbed = montecarlo_one_packet_loop (storage, &virt_packet, 1);
	  if ((virt_packet.nu < storage->spectrum_end_nu) &&
	      (virt_packet.nu > storage->spectrum_start_nu))
	    {
	      virt_id_nu =
		floor ((virt_packet.nu -
			storage->spectrum_start_nu) /
		       storage->spectrum_delta_nu);
	      storage->spectrum_virt_nu[virt_id_nu] +=
		virt_packet.energy * weight;
	    }
	}
    }
  return reabsorbed;
}

int64_t
montecarlo_one_packet_loop (storage_model_t * storage, rpacket_t * packet,
			    int64_t virtual_packet)
{
  rpacket_set_tau_event (packet, 0.0);
  rpacket_set_nu_line (packet, 0.0);
  rpacket_set_virtual_packet (packet, virtual_packet);
  rpacket_set_status (packet, TARDIS_PACKET_STATUS_IN_PROCESS);
  // Initializing tau_event if it's a real packet.
  if (virtual_packet == 0)
    {
      rpacket_reset_tau_event (packet);
    }
  // For a virtual packet tau_event is the sum of all the tau's that the packet passes.
  while (rpacket_get_status (packet) == TARDIS_PACKET_STATUS_IN_PROCESS)
    {
      // Check if we are at the end of line list.
      if (!rpacket_get_last_line (packet))
	{
	  rpacket_set_nu_line (packet,
			       storage->
			       line_list_nu[rpacket_get_next_line_id
					    (packet)]);
	}
      double distance;
      get_event_handler (packet, storage, &distance) (packet, storage,
						      distance);
      if (virtual_packet > 0 && rpacket_get_tau_event (packet) > 10.0)
	{
	  rpacket_set_tau_event (packet, 100.0);
	  rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
	}
    }
  if (virtual_packet > 0)
    {
      rpacket_set_energy (packet,
			  rpacket_get_energy (packet) * exp (-1.0 *
							     rpacket_get_tau_event
							     (packet)));
    }
  return rpacket_get_status (packet) ==
    TARDIS_PACKET_STATUS_REABSORBED ? 1 : 0;
}

tardis_error_t
rpacket_init (rpacket_t * packet, storage_model_t * storage, int packet_index,
	      int virtual_packet_flag)
{
  double nu_line;
  double current_r;
  double current_mu;
  double current_nu;
  double comov_current_nu;
  double current_energy;
  int64_t current_line_id;
  int current_shell_id;
  bool last_line;
  bool close_line;
  int recently_crossed_boundary;
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  current_nu = storage->packet_nus[packet_index];
  current_energy = storage->packet_energies[packet_index];
  current_mu = storage->packet_mus[packet_index];
  comov_current_nu = current_nu;
  current_shell_id = 0;
  current_r = storage->r_inner[0];
  current_nu =
    current_nu / (1 -
		  (current_mu * current_r * storage->inverse_time_explosion *
		   INVERSE_C));
  current_energy =
    current_energy / (1 -
		      (current_mu * current_r *
		       storage->inverse_time_explosion * INVERSE_C));
  if ((ret_val =
       line_search (storage->line_list_nu, comov_current_nu,
		    storage->no_of_lines,
		    &current_line_id)) != TARDIS_ERROR_OK)
    {
      return ret_val;
    }
  last_line = (current_line_id == storage->no_of_lines);
  recently_crossed_boundary = true;
  rpacket_set_nu (packet, current_nu);
  rpacket_set_mu (packet, current_mu);
  rpacket_set_energy (packet, current_energy);
  rpacket_set_r (packet, current_r);
  rpacket_set_current_shell_id (packet, current_shell_id);
  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  rpacket_set_close_line (packet, false);
  rpacket_set_recently_crossed_boundary (packet, recently_crossed_boundary);
  rpacket_set_virtual_packet_flag (packet, virtual_packet_flag);
  return ret_val;
}

/*
  Getter and setter methods.
*/

INLINE double
rpacket_get_nu (rpacket_t * packet)
{
  return packet->nu;
}

INLINE void
rpacket_set_nu (rpacket_t * packet, double nu)
{
  packet->nu = nu;
}

INLINE double
rpacket_get_mu (rpacket_t * packet)
{
  return packet->mu;
}

INLINE void
rpacket_set_mu (rpacket_t * packet, double mu)
{
  packet->mu = mu;
}

INLINE double
rpacket_get_energy (rpacket_t * packet)
{
  return packet->energy;
}

INLINE void
rpacket_set_energy (rpacket_t * packet, double energy)
{
  packet->energy = energy;
}

INLINE double
rpacket_get_r (rpacket_t * packet)
{
  return packet->r;
}

INLINE void
rpacket_set_r (rpacket_t * packet, double r)
{
  packet->r = r;
}

INLINE double
rpacket_get_tau_event (rpacket_t * packet)
{
  return packet->tau_event;
}

INLINE void
rpacket_set_tau_event (rpacket_t * packet, double tau_event)
{
  packet->tau_event = tau_event;
}

INLINE double
rpacket_get_nu_line (rpacket_t * packet)
{
  return packet->nu_line;
}

INLINE void
rpacket_set_nu_line (rpacket_t * packet, double nu_line)
{
  packet->nu_line = nu_line;
}

INLINE unsigned int
rpacket_get_current_shell_id (rpacket_t * packet)
{
  return packet->current_shell_id;
}

INLINE void
rpacket_set_current_shell_id (rpacket_t * packet,
			      unsigned int current_shell_id)
{
  packet->current_shell_id = current_shell_id;
}

INLINE unsigned int
rpacket_get_next_line_id (rpacket_t * packet)
{
  return packet->next_line_id;
}

INLINE void
rpacket_set_next_line_id (rpacket_t * packet, unsigned int next_line_id)
{
  packet->next_line_id = next_line_id;
}

INLINE bool
rpacket_get_last_line (rpacket_t * packet)
{
  return packet->last_line;
}

INLINE void
rpacket_set_last_line (rpacket_t * packet, bool last_line)
{
  packet->last_line = last_line;
}

INLINE bool
rpacket_get_close_line (rpacket_t * packet)
{
  return packet->close_line;
}

INLINE void
rpacket_set_close_line (rpacket_t * packet, bool close_line)
{
  packet->close_line = close_line;
}

INLINE int
rpacket_get_recently_crossed_boundary (rpacket_t * packet)
{
  return packet->recently_crossed_boundary;
}

INLINE void
rpacket_set_recently_crossed_boundary (rpacket_t * packet,
				       int recently_crossed_boundary)
{
  packet->recently_crossed_boundary = recently_crossed_boundary;
}

INLINE int
rpacket_get_virtual_packet_flag (rpacket_t * packet)
{
  return packet->virtual_packet_flag;
}

INLINE void
rpacket_set_virtual_packet_flag (rpacket_t * packet, int virtual_packet_flag)
{
  packet->virtual_packet_flag = virtual_packet_flag;
}

INLINE int
rpacket_get_virtual_packet (rpacket_t * packet)
{
  return packet->virtual_packet;
}

INLINE void
rpacket_set_virtual_packet (rpacket_t * packet, int virtual_packet)
{
  packet->virtual_packet = virtual_packet;
}

INLINE double
rpacket_get_d_boundary (rpacket_t * packet)
{
  return packet->d_boundary;
}

INLINE void
rpacket_set_d_boundary (rpacket_t * packet, double d_boundary)
{
  packet->d_boundary = d_boundary;
}

INLINE double
rpacket_get_d_electron (rpacket_t * packet)
{
  return packet->d_electron;
}

INLINE void
rpacket_set_d_electron (rpacket_t * packet, double d_electron)
{
  packet->d_electron = d_electron;
}

INLINE double
rpacket_get_d_line (rpacket_t * packet)
{
  return packet->d_line;
}

INLINE void
rpacket_set_d_line (rpacket_t * packet, double d_line)
{
  packet->d_line = d_line;
}

INLINE int
rpacket_get_next_shell_id (rpacket_t * packet)
{
  return packet->next_shell_id;
}

INLINE void
rpacket_set_next_shell_id (rpacket_t * packet, int next_shell_id)
{
  packet->next_shell_id = next_shell_id;
}

INLINE rpacket_status_t
rpacket_get_status (rpacket_t * packet)
{
  return packet->status;
}

INLINE void
rpacket_set_status (rpacket_t * packet, rpacket_status_t status)
{
  packet->status = status;
}

/* Other accessor methods. */

INLINE void
rpacket_reset_tau_event (rpacket_t * packet)
{
  rpacket_set_tau_event (packet, -log (rk_double (&mt_state)));
}
