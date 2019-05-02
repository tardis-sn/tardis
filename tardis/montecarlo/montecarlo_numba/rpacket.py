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
    ('close_line', boolean),
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
        self.current_shell_id = 0
        self.delta_shell_id = 0
        self.d_boundary = -1.0
        self.d_electron = -1.0
        self.d_line = 1.e99
        self.distance = 0.0
        self.nu_line = -1e99
        self.close_line = False
        self.last_line = False
        
        self.comov_nu_history = []
        self.radius_history = []
        self.move_dist_history = []
        self.next_interaction_history = []
        self.mu_history = []
        self.shell_id_history = []
        self.next_line_id_history = []
        
    def compute_distances(self, storage):
        """
        Compute all distances (d_line, d_boundary, ???), compare, 
        and set interaction
        
        Parameters
        ----------
        storage : [type]
            [description]
        """

        
        compute_distance2line(self, storage)
        compute_distance2boundary(self, storage)
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

    def move_packet(self, storage, distance):
        doppler_factor = self.get_doppler_factor(storage)
        r = self.r
        if (distance > 0.0):
            new_r = np.sqrt(r * r + distance * distance +
                             2.0 * r * distance * self.mu)
            self.mu = (self.mu * r + distance) / new_r
            self.r = new_r
        
        self.comov_nu_history.append(self.nu * doppler_factor)
        self.radius_history.append(self.r)
        self.move_dist_history.append(distance)
        self.next_interaction_history.append(self.next_interaction)
        self.mu_history.append(self.mu)
        self.shell_id_history.append(self.current_shell_id)
        self.next_line_id_history.append(self.next_line_id)
        

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
        self.next_line_id = next_line_id
        if self.next_line_id > (storage_model.no_of_lines - 1):
            self.last_line = True
        else:
            self.nu_line = storage_model.line_list_nu[self.next_line_id]
            self.last_line = False
        
        ##### FIXME Add close line initializer in a sensible  - think about this!1
        #self.set_close_line(storage_model)

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
        if self.distance == 0.0:
            self.set_close_line(storage_model)
        next_line_id = self.next_line_id
        storage_model.line_lists_tau_sobolevs
        tau_line = storage_model.line_lists_tau_sobolevs[next_line_id, 
                        self.current_shell_id]
        tau_continuum = 0.0
        tau_combined = tau_line + tau_continuum
        
        if (next_line_id + 1) == storage_model.no_of_lines:
            self.last_line = True
        if (self.tau_event < tau_combined): # Line absorption occurs
            self.move_packet(storage_model, self.distance)
            self.transform_energy(storage_model)
            self.line_emission(storage_model)
        else:
            self.tau_event -= tau_line
            self.next_line_id = next_line_id + 1
            if not self.last_line:
                self.nu_line = storage_model.line_list_nu[self.next_line_id]
        


    def line_emission(self, storage_model):
        emission_line_id = self.next_line_id
        inverse_doppler_factor = 1 / self.get_doppler_factor(storage_model)
        self.nu = storage_model.line_list_nu[emission_line_id] * inverse_doppler_factor 
        
        self.next_line_id = emission_line_id + 1
        self.nu_line = storage_model.line_list_nu[self.next_line_id]
        self.tau_event = get_tau_event()
        self.set_close_line(storage_model)


    def set_close_line(self, storage_model, line_diff_threshold=1e-7):
        """
        The Packet currently sits on next_line_id - 1 and we are checking if the
        next line is closer than 1e-7 to set the close_line attribute

        Parameters
        ----------
        storage_model : StorageModel
        """
        return
        if self.last_line:
            self.close_line = False
            return
        
        frac_nu_diff = ((storage_model.line_list_nu[self.next_line_id] -
                        storage_model.line_list_nu[self.next_line_id - 1]) /
                         storage_model.line_list_nu[self.next_line_id])
        
        if not self.last_line and frac_nu_diff < line_diff_threshold:
            self.close_line = True
        else:
            self.close_line = False
#{
#  if (!rpacket_get_last_line (packet) &&
#      fabs (storage->line_list_nu[rpacket_get_next_line_id (packet)] -
#            rpacket_get_nu_line (packet)) < (rpacket_get_nu_line (packet)*
#                                             1e-7))
#    {
#      rpacket_set_close_line (packet, true);
#    }
#}
