#include "rpacket.h"
#include "storage.h"

tardis_error_t
rpacket_init (rpacket_t * packet, storage_model_t * storage, int packet_index,
              int virtual_packet_flag)
{
  int64_t current_line_id;
  tardis_error_t ret_val = TARDIS_ERROR_OK;
  double current_nu = storage->packet_nus[packet_index];
  double current_energy = storage->packet_energies[packet_index];
  double current_mu = storage->packet_mus[packet_index];
  double comov_current_nu = current_nu;
  int current_shell_id = 0;
  double current_r = storage->r_inner[0];
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
  bool last_line = (current_line_id == storage->no_of_lines);
  rpacket_set_nu (packet, current_nu);
  rpacket_set_mu (packet, current_mu);
  rpacket_set_energy (packet, current_energy);
  rpacket_set_r (packet, current_r);
  rpacket_set_current_shell_id (packet, current_shell_id);
  rpacket_set_next_line_id (packet, current_line_id);
  rpacket_set_last_line (packet, last_line);
  rpacket_set_close_line (packet, false);
  rpacket_set_recently_crossed_boundary (packet, true);
  rpacket_set_virtual_packet_flag (packet, virtual_packet_flag);
  return ret_val;
}
