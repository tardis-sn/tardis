import numpy as np
from astropy import units as u, constants as const
class Model():
    pass

class HomologousRadial1D(Model):
    def __init__(self, v_outer):
        self.v_inner = np.zeros_like(v_outer)
        self.v_inner[1:] = v_outer[:-1]