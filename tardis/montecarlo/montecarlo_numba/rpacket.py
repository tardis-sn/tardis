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
    ('tau_event', float64),
    ('nu_line', float64),
    ('last_line', boolean),
    ('next_line_id', int64),
    ('d_line', float64),  # Distance to line event. 
    ('d_electron', float64), #/**< Distance to line event. */
    ('d_boundary', float64), # distance to boundary 
    ('current_shell_id', int64),
    ('delta_shell_id', int64),
    ('next_interaction', int64),
    ('status', int64),
    ('distance', float64)
]

@njit(**njit_dict)
def get_tau_event():
    return  np.random.exponential()


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
        self.tau_event = get_tau_event()
        #self.nu_line = line_search() # TODO: Implement this
        self.current_shell_id = 0
        self.delta_shell_id = 0
        self.d_boundary = -1.0
        self.d_electron = -1.0
        self.d_line = 1.e99
        self.distance = 0.0
        self.nu_line = -1e99
    
    def compute_distances(self, storage):
        compute_distance2line(self, storage)
        compute_distance2boundary(self, storage)
        print('d_boundary', self.d_boundary, 'd_line', self.d_line)
        if self.d_boundary < self.d_line:
            next_interaction = BOUNDARY
            self.distance = self.d_boundary
        else:
            next_interaction = LINE
            self.distance = self.d_line

        self.next_interaction = next_interaction

    
    def get_doppler_factor(self, storage):
        beta = self.r * storage.inverse_time_explosion / C_SPEED_OF_LIGHT
        
        return 1.0 - self.mu * beta
        #else:
        #    return (1.0 - self.mu * beta) / np.sqrt(1 - beta*beta)

    def move_packet(self, storage, distance):
        doppler_factor = self.get_doppler_factor(storage)
        r = self.r
        if (distance > 0.0):
            new_r = np.sqrt(r * r + distance * distance +
                             2.0 * r * distance * self.mu)
            self.mu = (self.mu * r + distance) / new_r
            self.r = new_r

    def move_packet_across_shell_boundary(self, storage):
        self.move_packet(storage, self.distance)

        get_tau_event()
        if ((self.current_shell_id < storage.no_of_shells - 1 and self.delta_shell_id == 1) 
            or (self.current_shell_id > 0 and self.delta_shell_id == -1)):
            self.current_shell_id += self.delta_shell_id
        elif self.delta_shell_id == 1:
            self.status = EMITTED
        else:
            self.status = REABSORBED
    
    def set_line(self, storage_model):
        inverse_line_list_nu = storage_model.line_list_nu[::-1]
        doppler_factor = self.get_doppler_factor(storage_model)
        comov_nu = self.nu * doppler_factor
        next_line_id = storage_model.no_of_lines - np.searchsorted(inverse_line_list_nu, comov_nu)
        print('in set line comov nu', comov_nu, 'next_line_nu', storage_model.line_list_nu[next_line_id-1:next_line_id+2])
        #print('packet nu', self.nu, 'next_line_id', next_line_id, 'next_line_nu', storage_model.line_list_nu[next_line_id-1:next_line_id+2])
        #self.next_line_id = 3000
        self.next_line_id = next_line_id
        print('in set_line nextid and len', self.next_line_id, len(storage_model.line_list_nu))
        if self.next_line_id > (storage_model.no_of_lines - 1):
            self.last_line = True
        else:
            self.nu_line = storage_model.line_list_nu[self.next_line_id]
            self.last_line = False

 # if (rpacket_get_virtual_packet (packet) > 0)
 #   {
#      rpacket_set_tau_event (packet,
#                             rpacket_get_tau_event (packet) + tau_line);
#      rpacket_set_next_line_id (packet, next_line_id + 1);
#      test_for_close_line (packet, storage);
#    } 
    def transform_energy(self, storage_model):
        """
        Transform from the LabFrame to the ComovingFrame. Then change the angle 
        and transform back conserving energy in the ComovingFrame.        
        """
        old_doppler_factor = self.get_doppler_factor(storage_model)
        self.mu = get_random_mu()
        inverse_doppler_factor = 1. / self.get_doppler_factor(storage_model)
        comov_energy = self.energy * old_doppler_factor
        self.energy = comov_energy * inverse_doppler_factor

    def line_scatter(self, storage_model):
        #print('Line Scattering')
        next_line_id = self.next_line_id
        #tau_line = storage_model.line_lists_tau_sobolevs[next_line_id, self.current_shell_id]
        tau_line = 3.0
        # TODO: Fixme
        tau_continuum = 0.0
        tau_combined = tau_line + tau_continuum
        
        if (next_line_id + 1) == storage_model.no_of_lines:
            self.last_line = True
        if (self.tau_event < tau_combined): # Line absorption occurs
            self.move_packet(storage_model, self.distance)
            self.transform_energy(storage_model)
            self.line_emission(storage_model)
            #print('rpacket scattered at', self.nu)
        else:
            self.tau_event -= tau_line
            self.next_line_id = next_line_id + 1
            # ???
            self.nu_line = storage_model.line_list_nu[self.next_line_id] 
            #test_for_close_line (packet, storage);

    def line_emission(self, storage_model):
        emission_line_id = self.next_line_id
        inverse_doppler_factor = 1 / self.get_doppler_factor(storage_model)
        self.nu = storage_model.line_list_nu[emission_line_id] * inverse_doppler_factor 
        #self.nu_line =  storage->line_list_nu[emission_line_id]);
        self.next_line_id = emission_line_id + 1
        self.tau_event = get_tau_event()
#void test_for_close_line (rpacket_t * packet, const storage_model_t * storage)
#{
#  if (!rpacket_get_last_line (packet) &&
#      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
#            rpacket_get_nu_line (packet)) < (rpacket_get_nu_line (packet)*
#                                             1e-7))
#    {
#      rpacket_set_close_line (packet, true);
#    }
#}
