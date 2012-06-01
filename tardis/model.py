#File contains model class

import constants
import numpy as np
import os
import plasma
miss_distance = 1e99

w7model_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'w7model.dat'))

class Model(object):
    pass


class MultiZoneRadial(Model):
    @classmethod
    def from_w7(cls, current_time):
        velocities, densities = np.loadtxt(w7model_file, usecols=(1, 2), unpack=1)
        velocities *= 1e5
        densities = 10**densities
        time_0 = 0.000231481 * 86400
        return cls(velocities[100::5], densities[100::5], current_time, time_0)
        
    
    def __init__(self, velocities, densities, current_time, time_0):
        #constructing zones
        self.v_inner = velocities[:-1]
        self.v_outer = velocities[1:]
        
        self.time_0 = time_0
        self.current_time = current_time
        
        self.r_inner = self.v_inner * current_time
        self.r_outer = self.v_outer * current_time
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)
        self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0]**2/ self.r_middle**2))
        self.densities_middle = (current_time / time_0)**-3 * densities[1:]
    
    def set_line_list(self, line_list):
        self.line_list = line_list
    
    def set_atomic_model(self, atomic_model):
        self.atomic_model = atomic_model
    
    def read_abundances_uniform(self, abundances):
        self.abundances = np.empty((self.r_inner.size, abundances.size))
        for i in range(self.r_inner.size):
            self.abundances[i] = abundances.copy()
    
    def _initialize_temperature(self, t_rad):
        self.t_rads = np.ones(self.r_inner.size) * t_rad
    
    def initialize_plasmas(self, t_rad):
        self._initialize_temperature(t_rad)
        self.electron_densities = np.empty(self.r_inner.size, dtype=np.float64)
        self.plasmas = []
        for i, (current_abundance, current_t, current_w, current_density) in enumerate(zip(self.abundances, self.t_rads, self.ws, self.densities_middle)):
            current_plasma = plasma.NebularPlasma.from_model(current_abundance, current_density, self.atomic_model, abundance_mode='array')
            current_plasma.update_radiationfield(t_rad=current_t, w=current_w)
            self.plasmas.append(
                current_plasma
                )
            self.electron_densities[i] = current_plasma.electron_density
    
    
    def calculate_tau_sobolevs(self):
        tau_sobolevs = np.empty((self.r_inner.size, len(self.line_list)))
        for i, plasma in enumerate(self.plasmas):
            tau_sobolevs[i] = plasma.calculate_tau_sobolev(self.line_list, self.current_time)
        return tau_sobolevs
    
    def calculate_transition_probabilities(self, tau_sobolevs):
        transition_probabilities = np.empty((self.r_inner.size, self.atomic_model.macro_atom.target_level_total.size))
        for i, (t_rad, w, tau_sobolev) in enumerate(zip(self.t_rads, self.ws, tau_sobolevs)):
            transition_probabilities[i] = self.atomic_model.macro_atom.calculate_transition_probabilities(tau_sobolev, t_rad, w)
        return transition_probabilities
    
    def update_model(self, t_rads, ws):
        self.t_rads = t_rads
        self.ws = ws
        for i, (w, t_rad, current_plasma) in enumerate(zip(self.ws, self.t_rads, self.plasmas)):
            current_plasma.update_radiationfield(t_rad, w)
            self.electron_densities[i] = current_plasma.electron_density
            
