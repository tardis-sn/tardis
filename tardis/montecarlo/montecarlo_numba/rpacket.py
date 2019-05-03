import numpy as np
from enum import Enum
from numba import int64, float64, boolean
from numba import jitclass, njit
from tardis.montecarlo.montecarlo_numba.compute_distance import (compute_distance2boundary,
                                                                 compute_distance2continuum,
                                                                 compute_distance2line)
from tardis.montecarlo.montecarlo_numba import njit_dict
from astropy import constants as const

C_SPEED_OF_LIGHT = const.c.to('cm/s').value
MISS_DISTANCE = 1e99

#class PacketStatus(Enum):
IN_PROCESS = 0
EMITTED = 1
REABSORBED = 2

#class InteractionType(Enum):
ESCATTERING = 0
BOUNDARY = 1
LINE = 2

rpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64),
    ('next_line_id', int64),
    ('current_shell_id', int64),
    ('status', int64),
]

@njit(**njit_dict)
def calculate_distance_boundary(r, mu, r_inner, r_outer):
    delta_shell = 0
    if (mu > 0.0):
        # direction outward
        distance = np.sqrt(r_outer**2 + ((mu**2 - 1.0) * r**2)) - (r * mu)
        delta_shell = 1
    else:
        # going inward
        check = r_inner**2 + (r**2 * (mu**2 - 1.0))

        if (check >= 0.0):
            # hit inner boundary 
            distance = -r * mu - np.sqrt(check)
            delta_shell = -1
        else:
            # miss inner boundary 
            distance = np.sqrt(r_outer**2 + ((mu**2 - 1.0) * r**2)) - (r * mu)
            delta_shell = 1
    
    return distance, delta_shell

@njit(**njit_dict)
def calculate_distance_line(nu, comov_nu, nu_line, ct):
    if nu_line == 0.0:
        return MISS_DISTANCE
    nu_diff = comov_nu - nu_line
    if np.abs(nu_diff / comov_nu) < 1e-7:
            nu_diff = 0.0
    if nu_diff >= 0:                    
        return (nu_diff / nu) * ct
    else:
        #return np.abs((nu_diff/self.nu) * ct)
        raise Exception

@njit(**njit_dict)
def get_doppler_factor(r, mu, inverse_time_explosion):
    beta = (r * inverse_time_explosion) / C_SPEED_OF_LIGHT
    return 1.0 - mu * beta


@njit(**njit_dict)
def get_random_mu():
    return 2.0 * np.random.random() - 1.0


@jitclass(rpacket_spec)
class RPacket(object):
    def __init__(self, r, mu, nu, energy):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = IN_PROCESS
        self.next_line_id = -1
        
    def calculate_distance_continuum(self, storage):    
        packet.d_electron = storage.inverse_electron_densities[packet.current_shell_id] * \
                storage.inverse_sigma_thomson * packet.tau_event
        
    def trace_packet(self, storage_model):
        r_inner = storage_model.r_inner[self.current_shell_id]
        r_outer = storage_model.r_outer[self.current_shell_id]
        
        distance = 0.0

        distance_boundary, delta_shell = calculate_distance_boundary(self.r, self.mu, r_inner, r_outer)
        
        #defining start for stuff
        cur_line_id = self.next_line_id
        nu_line = 0.0
        #defining taus
        tau_event = np.random.exponential()
        tau_trace_line = 0.0
        tau_trace_line_combined = 0.0
        #Calculating doppler factor
        doppler_factor = get_doppler_factor(self.r, self.mu, storage_model.inverse_time_explosion)
        comov_nu = self.nu * doppler_factor
        distance_trace = 0.0
        d_continuum = MISS_DISTANCE
        
        while True:
            if cur_line_id < storage_model.no_of_lines: # last_line
                nu_line = storage_model.line_list_nu[cur_line_id]
                tau_trace_line += storage_model.line_lists_tau_sobolevs[cur_line_id, 
                        self.current_shell_id]
            else:
                nu_line = 0.0
                tau_trace_line = 0.0
                interaction_type = BOUNDARY  # FIXME: does not work for e-scattering
                distance = distance_boundary
                break
            
            tau_trace_line_combined += tau_trace_line
            distance_trace = calculate_distance_line(self.nu, comov_nu, nu_line, storage_model.ct)
            tau_trace_combined = tau_trace_line_combined + 0 #tau_trace_electron electron scattering
            
            if distance_trace > distance_boundary:
                interaction_type = BOUNDARY # BOUNDARY
                self.next_line_id = cur_line_id
                distance = distance_boundary
                break
            
            if distance_trace > d_continuum:
                interaction_type = 10 #continuum
                #break
            if tau_trace_combined > tau_event:
                interaction_type = LINE #Line
                self.next_line_id = cur_line_id
                distance = distance_trace
                break
            
            cur_line_id += 1
            
            #loop_checker +=1
            #if loop_checker > 10000:
            #    raise Exception
         
        return distance, interaction_type, delta_shell

    def move_packet(self, distance):
        """Move packet a distance and recalculate the new angle mu
        
        Parameters
        ----------
        distance : float
            distance in cm
        """

        r = self.r
        if (distance > 0.0):
            new_r = np.sqrt(r**2 + distance**2 +
                             2.0 * r * distance * self.mu)
            self.mu = (self.mu * r + distance) / new_r
            self.r = new_r

    def move_packet_across_shell_boundary(self, distance, delta_shell, no_of_shells):
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

        self.move_packet(distance)
        if ((self.current_shell_id < no_of_shells - 1 and delta_shell == 1) 
            or (self.current_shell_id > 0 and delta_shell == -1)):
            self.current_shell_id += delta_shell
        elif delta_shell == 1:
            self.status = EMITTED
        else:
            self.status = REABSORBED
    
    def initialize_line_id(self, storage_model):
        inverse_line_list_nu = storage_model.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(self.r, self.mu, storage_model.inverse_time_explosion)
        comov_nu = self.nu * doppler_factor
        next_line_id = storage_model.no_of_lines - np.searchsorted(inverse_line_list_nu, comov_nu)
        self.next_line_id = next_line_id

    def transform_energy(self, storage_model):
        """
        Transform from the LabFrame to the ComovingFrame. Then change the angle 
        and transform back conserving energy in the ComovingFrame.        
        """
        old_doppler_factor = get_doppler_factor(self.r, self.mu, storage_model.inverse_time_explosion)
        self.mu = get_random_mu()
        inverse_new_doppler_factor = 1. / get_doppler_factor(self.r, self.mu, storage_model.inverse_time_explosion)
        comov_energy = self.energy * old_doppler_factor
        self.energy = comov_energy * inverse_new_doppler_factor

    def line_emission(self, storage_model):
        emission_line_id = self.next_line_id
        inverse_doppler_factor = 1 / get_doppler_factor(self.r, self.mu, storage_model.inverse_time_explosion)
        self.nu = storage_model.line_list_nu[emission_line_id] * inverse_doppler_factor 
        self.next_line_id = emission_line_id + 1