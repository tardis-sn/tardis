#include "rpacket.h"
#include "storage.h"

rk_state mt_state;

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

inline int rpacket_get_id (rpacket_t * packet)
{
  return packet->id;
}

inline void rpacket_set_id (rpacket_t * packet, int id)
{
  packet->id = id;
}

/* Other accessor methods. */

INLINE void
rpacket_reset_tau_event (rpacket_t * packet)
{
  rpacket_set_tau_event (packet, -log (rk_double (&mt_state)));
}
