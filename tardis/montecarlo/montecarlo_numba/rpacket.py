import numpy as np
from enum import Enum
from numba import int64, float64
from numba import jitclass, njit
from tardis.montecarlo.montecarlo_numba.compute_distance import (compute_distance2boundary,
                                                                 compute_distance2continuum)

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
 #   ('d_line', float64)  # Distance to electron event. 
    ('d_electron', float64), #/**< Distance to line event. */
    ('d_boundary', float64), # distance to boundary 
    ('current_shell_id', int64),
    ('delta_shell_id', int64),
    ('next_interaction', int64),
    ('status', int64),
    ('distance', float64)
]

@njit
def get_tau_event():
    return  np.random.exponential()


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
        self.distance = 0.0
        #self.nu_line = nu_line
    
    def compute_distances(self, storage):
        compute_distance2continuum(self, storage)
        compute_distance2boundary(self, storage)
        if self.d_boundary < self.d_electron:
            next_interaction = BOUNDARY
            self.distance = self.d_boundary
        else:
            next_interaction = ESCATTERING
            # TODO: Fixme
            self.distance = self.d_boundary

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
