from numba import float64, int64
import numpy as np

from
from tardis.montecarlo.montecarlo_numba.rpacket import (
    calculate_distance_boundary, calculate_distance_electron,
    get_doppler_factor, calculate_distance_line, calculate_tau_electron)

vpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64),
    ('next_line_id', int64),
    ('current_shell_id', int64),
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
        cur_inverse_electron_density = 1 / cur_electron_density
        distance_electron = calculate_distance_electron(
            cur_inverse_electron_density, tau_event)


        #Calculating doppler factor
        doppler_factor = get_doppler_factor(self.r, self.mu, 
                                        storage_model.inverse_time_explosion)
        comov_nu = self.nu * doppler_factor
        distance_trace = 0.0
        last_line = False

        tau_trace_electron = calculate_tau_electron(cur_electron_density,
                                                    distance_trace)

        tau_trace_combined = tau_trace_electron

        while self.status == INPROCESS:
            if cur_line_id < storage_model.no_of_lines:  # not last_line
                nu_line = storage_model.line_list_nu[cur_line_id]
                tau_trace_line = storage_model.line_lists_tau_sobolevs[
                    cur_line_id, self.current_shell_id]
            else:
                last_line = True
                self.next_line_id = cur_line_id
                break
            
            tau_trace_combined += tau_trace_line
            distance_trace = calculate_distance_line(self.nu, comov_nu, nu_line, 
                                                        storage_model.ct)

            tau_trace_combined = tau_trace_line_combined + tau_trace_electron

            if (distance_boundary <= distance_trace):
                current_shell_id += delta
                self.move_packet_across_shell_boundary(distance_boundary)
                distance_trace =
                distance_boundary, delta_shell = calculate_distance_boundary(
                    self.r, self.mu, r_inner, r_outer)

                ## distance_boundary
                #unless shell 
                interaction_type = InteractionType.BOUNDARY # BOUNDARY
                self.next_line_id = cur_line_id
                distance = distance_boundary
                break
            
            if (distance_electron < distance_trace) and (distance_electron < distance_boundary):
                interaction_type = InteractionType.ESCATTERING
                distance = distance_electron
                self.next_line_id = cur_line_id
                break
                        
            cur_line_id += 1
        if not last_line:            
            return distance, interaction_type, delta_shell
        else:
            if distance_electron < distance_boundary:
                return distance_electron, InteractionType.ESCATTERING, delta_shell
            else:
                return distance_boundary, InteractionType.BOUNDARY, delta_shell

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
