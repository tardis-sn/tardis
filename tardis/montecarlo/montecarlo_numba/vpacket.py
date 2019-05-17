from numba import float64, int64, jitclass
import numpy as np

from
from tardis.montecarlo.montecarlo_numba.rpacket import (
    calculate_distance_boundary, get_doppler_factor, calculate_distance_line,
    calculate_tau_electron, PacketStatus)

vpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64),
    ('next_line_id', int64),
    ('current_shell_id', int64),
    ('status', int64)
]

@jitclass(vpacket_spec)
class VPacket(object):
    def __init__(self, r, mu, nu, energy):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = -1
        self.next_line_id = -1
        self.status = -1

    def trace_vpacket(self, storage_model):
        
        r_inner = storage_model.r_inner[self.current_shell_id]
        r_outer = storage_model.r_outer[self.current_shell_id]
        
        distance = 0.0

        distance_boundary, delta_shell = calculate_distance_boundary(
            self.r, self.mu, r_inner, r_outer)
        
        #defining start for line interaction
        cur_line_id = self.next_line_id
        nu_line = 0.0

        #defining taus
        tau_event = np.random.exponential()
        tau_trace_line = 0.0
        tau_trace_line_combined = 0.0
        
        #e scattering initialization

        cur_electron_density = storage_model.electron_densities[
            self.current_shell_id]

        #Calculating doppler factor
        doppler_factor = get_doppler_factor(self.r, self.mu, 
                                        storage_model.inverse_time_explosion)
        last_line = False

        tau_trace_electron = calculate_tau_electron(cur_electron_density,
                                                    distance_boundary)

        tau_trace_combined = tau_trace_electron

        while True:
            if cur_line_id < storage_model.no_of_lines:  # not last_line
                nu_line = storage_model.line_list_nu[cur_line_id]
                tau_trace_line = storage_model.line_lists_tau_sobolevs[
                    cur_line_id, self.current_shell_id]
                tau_trace_combined += tau_trace_line
                distance_trace = calculate_distance_line(self.nu, comov_nu,
                                                         nu_line,
                                                         storage_model.ct)
                cur_line_id += 1
            else:
                distance_trace = distance_boundary + 1.0
            



            if distance_boundary <= distance_trace:
                self.current_shell_id += delta_shell
                self.move_packet_across_shell_boundary(distance_boundary)
                if not self.status == PacketStatus.IN_PROCESS:
                    break

                comov_nu = self.nu * doppler_factor
                tau_electron = calculate_tau_electron(cur_electron_density,
                                                      distance_boundary)
                tau_trace_combined += tau_electron

            return tau_trace_combined

    def move_vpacket(self, distance):
        """Move packet a distance and recalculate the new angle mu

        Parameters
        ----------
        distance : float
            distance in cm
        """
        r = self.r
        if (distance > 0.0):
            new_r = np.sqrt(r ** 2 + distance ** 2 +
                            2.0 * r * distance * self.mu)
            self.mu = (self.mu * r + distance) / new_r
            self.r = new_r

    def move_packet_across_shell_boundary(self, distance, delta_shell,
                                          no_of_shells):
        """
        Move packet across shell boundary - realizing if we are still in the simulation or have
        moved out through the inner boundary or outer boundary and updating packet
        status.

        Parameters
        ----------
        distance : float
            distance to move to shell boundary

        delta_shell: int
            is +1 if moving outward or -1 if moving inward

        no_of_shells: int
            number of shells in TARDIS simulation
        """

        if ((self.current_shell_id < no_of_shells - 1 and delta_shell == 1)
                or (self.current_shell_id > 0 and delta_shell == -1)):
            self.current_shell_id += delta_shell
        elif delta_shell == 1:
            self.status = PacketStatus.EMITTED
        else:
            self.status = PacketStatus.REABSORBED


def create_single_vpacket():

def create_volley_vpacket():

    if cur_r_packet.r > r_innerst:
        mu_min = -np.sqrt(1 - r_innerst / cur_r_packet.r) ** 2
    else:
        mu_min = 0.0

    double
    mu_bin = (1.0 - mu_min) / rpacket_get_virtual_packet_flag(packet);
    rpacket_set_mu( & virt_packet, mu_min + (i + rk_double(mt_state)) * mu_bin);
    {
      if ((rpacket_get_nu (packet) > storage->spectrum_virt_start_nu) && (rpacket_get_nu(packet) < storage->spectrum_virt_end_nu))
        {
          for (int64_t i = 0; i < rpacket_get_virtual_packet_flag (packet); i++)
            {
              double weight;
              rpacket_t virt_packet = *packet;
              double mu_min;
              if (rpacket_get_r(&virt_packet) > storage->r_inner[0])
                {
                  mu_min =
                    -1.0 * sqrt (1.0 -
                                 (storage->r_inner[0] / rpacket_get_r(&virt_packet)) *
                                 (storage->r_inner[0] / rpacket_get_r(&virt_packet)));

                  if (storage->full_relativity)
                    {
                      // Need to transform the angular size of the photosphere into the CMF
                      mu_min = angle_aberration_LF_to_CMF (&virt_packet, storage, mu_min);
                    }
                }
              else
                {
                  mu_min = 0.0;
                }
              double mu_bin = (1.0 - mu_min) / rpacket_get_virtual_packet_flag (packet);
              rpacket_set_mu(&virt_packet,mu_min + (i + rk_double (mt_state)) * mu_bin);
              switch (virtual_mode)
                {
                case -2:
                  weight = 1.0 / rpacket_get_virtual_packet_flag (packet);
                  break;
                case -1:
                  weight =
                    2.0 * rpacket_get_mu(&virt_packet) /
                    rpacket_get_virtual_packet_flag (packet);
                  break;
                case 1:
                  weight =
                    (1.0 -
                     mu_min) / 2.0 / rpacket_get_virtual_packet_flag (packet);
                  break;
                default:
                  fprintf (stderr, "Something has gone horribly wrong!\n");
                  // FIXME MR: we need to somehow signal an error here
                  // I'm adding an exit() here to inform the compiler about the impossible path
                  exit(1);
                }
              angle_aberration_CMF_to_LF (&virt_packet, storage);
              double doppler_factor_ratio =
                rpacket_doppler_factor (packet, storage) /
                rpacket_doppler_factor (&virt_packet, storage);
              rpacket_set_energy(&virt_packet,
                                 rpacket_get_energy (packet) * doppler_factor_ratio);
              rpacket_set_nu(&virt_packet,rpacket_get_nu (packet) * doppler_factor_ratio);
              reabsorbed = montecarlo_one_packet_loop (storage, &virt_packet, 1, mt_state);
#ifdef WITH_VPACKET_LOGGING
#ifdef WITHOPENMP
#pragma omp critical
                {
#endif // WITHOPENMP
                  if (storage->virt_packet_count >= storage->virt_array_size)
                    {
                      storage->virt_array_size *= 2;
                      storage->virt_packet_nus = safe_realloc(storage->virt_packet_nus, sizeof(double) * storage->virt_array_size);
                      storage->virt_packet_energies = safe_realloc(storage->virt_packet_energies, sizeof(double) * storage->virt_array_size);
                      storage->virt_packet_last_interaction_in_nu = safe_realloc(storage->virt_packet_last_interaction_in_nu, sizeof(double) * storage->virt_array_size);
                      storage->virt_packet_last_interaction_type = safe_realloc(storage->virt_packet_last_interaction_type, sizeof(int64_t) * storage->virt_array_size);
                      storage->virt_packet_last_line_interaction_in_id = safe_realloc(storage->virt_packet_last_line_interaction_in_id, sizeof(int64_t) * storage->virt_array_size);
                      storage->virt_packet_last_line_interaction_out_id = safe_realloc(storage->virt_packet_last_line_interaction_out_id, sizeof(int64_t) * storage->virt_array_size);
                    }
                  storage->virt_packet_nus[storage->virt_packet_count] = rpacket_get_nu(&virt_packet);
                  storage->virt_packet_energies[storage->virt_packet_count] = rpacket_get_energy(&virt_packet) * weight;
                  storage->virt_packet_last_interaction_in_nu[storage->virt_packet_count] = storage->last_interaction_in_nu[rpacket_get_id (packet)];
                  storage->virt_packet_last_interaction_type[storage->virt_packet_count] = storage->last_interaction_type[rpacket_get_id (packet)];
                  storage->virt_packet_last_line_interaction_in_id[storage->virt_packet_count] = storage->last_line_interaction_in_id[rpacket_get_id (packet)];
                  storage->virt_packet_last_line_interaction_out_id[storage->virt_packet_count] = storage->last_line_interaction_out_id[rpacket_get_id (packet)];
                  storage->virt_packet_count += 1;
#ifdef WITHOPENMP
                }
#endif // WITHOPENMP
#endif // WITH_VPACKET_LOGGING
              if ((rpacket_get_nu(&virt_packet) < storage->spectrum_end_nu) &&
                  (rpacket_get_nu(&virt_packet) > storage->spectrum_start_nu))
                {
#ifdef WITHOPENMP
#pragma omp critical
                    {
#endif // WITHOPENMP
                      int64_t virt_id_nu =
                        floor ((rpacket_get_nu(&virt_packet) -
                                storage->spectrum_start_nu) /
                               storage->spectrum_delta_nu);
                      storage->spectrum_virt_nu[virt_id_nu] +=
                        rpacket_get_energy(&virt_packet) * weight;
#ifdef WITHOPENMP
                    }
#endif // WITHOPENMP
                }
            }
        }
      else
        {
          return 1;
        }