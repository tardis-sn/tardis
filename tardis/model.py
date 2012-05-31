#File contains model class

import constants
import numpy as np
import os
miss_distance = 1e99

w7model_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'w7model.dat'))

class Model(object):
    pass


class MultiZoneRadial(Model):
    @classmethod
    def from_w7(cls, current_time):
        velocities, densities = loadtxt(w7model_file, usecols=(1, 2), unpack=1)
        velocities *= 1e5
        densities = 10**densities
        time_0 = 0.000231481 * 86400
        cls(velocities[100::5], densities[100::5], current_time, time_0)
        
    
    def __init__(self, velocities, densities, current_time, time_0):
        #constructing zones
        self.v_inner = velocities[:-1]
        self.v_outer = velocities[1:]
        
        
        self.r_inner = self.v_inner * time_exp
        self.r_outer = self.v_outer * time_exp
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)
        self.initial_w = 0.5 * (1 - np.sqrt(1 - self.r_inner[0]**2/ self.r_middle**2))
        self.densities_middle = (current_time / time_0)**-3 * self.densities[1:]
        
    def read_abundances_uniform(self, abundances):
        self.abundances = np.empty((self.r_inner.size, self.abundances.size))
        for line in abundances:
            self.abundances = line
            