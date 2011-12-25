#File contains model class

import constants
import numpy as np

miss_distance = 1e99

class Model(object):
    pass



class OneZoneBaseModel(Model):
    def __init__(self, r_inner, r_outer, v_outer):
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
        self.homol_coeff = self.v_outer / self.r_outer
    
    def compute_distance2outer(self, r, mu):
        #compute distance to the outer layer
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
        cur_v = r * self.homol_coeff
        #beta = np.abs(mu) * cur_v / constants.c
        nu_cmf = nu * (1. - (mu * r * self.homol_coeff / constants.c))
        #nu_cmf = nu * np.sqrt((1-beta)/(1+beta))
        
        #this_nucmf = this_freq * (1. - (this_mu * v1 * this_r / r1 / CLIGHT));
        if nu_cmf < nu_line:
            return miss_distance
        else:
            return ((nu_cmf - nu_line) / nu) * constants.c / self.homol_coeff        
    