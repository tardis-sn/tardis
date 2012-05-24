#File contains model class

import constants
import numpy as np
import c_distance
miss_distance = 1e99

class Model(object):
    pass



class OneZoneBaseModel(Model):
    def __init__(self, r_inner, r_outer, v_outer, ne):
        """Initialize the model
        Parameters
        r_inner : float
            inner radius
        r_outer : float
            outer radius
        v_inner : float
            velocity at inner radius
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_outer = v_outer
        self.t_exp = self.r_outer / self.v_outer
        self.inverse_t_exp = 1 / self.t_exp
        self.v_inner = self.r_inner / self.t_exp
        self.ne = ne
        self. inverse_ne = 1 / self.ne
        
    def compute_distance2outer(self, r, mu):
        #compute distance to the outer layer
        return c_distance.distance2outer(self.r_outer, r, mu)
        d = np.sqrt(self.r_outer**2 + ((mu**2 - 1.) * r**2)) - (r * mu)
        return d
    
    def compute_distance2inner(self, r, mu):
        #compute distance to the inner layer
        #check if intersection is possible?
        check = self.r_inner**2 + (r**2 * (mu**2 - 1.))
        if check < 0:
            return miss_distance
        else:
            if mu < 0:
               return -r * mu - np.sqrt(check)
            else:
                return miss_distance
    
    def compute_distance2line(self, r, mu, nu, nu_line):
        #computing distance to line

        nu_cmf = nu * (1. - (mu * r * self.inverse_t_exp * constants.inverse_c))

        if nu_cmf < nu_line:
            return miss_distance
        else:
            return ((nu_cmf - nu_line) / nu) * constants.c * self.t_exp
    
    def compute_distance2electron(self, r, mu, tau_event):
        return tau_event * self.inverse_ne * constants.inverse_sigma_thomson

class OneDMultiZoneModel(Model):
    pass
    