#creating photons

import numpy
import constants
class PhotonSource(object):
    pass

class SimplePhotonSource(PhotonSource):
    #Draw uniform randomly simple photon packets
    
    @classmethod
    def from_wavelength(cls, wl_start, wl_end, seed=250819801106):
        """Initializing from wavelength
        Parameters
        ----------
        
        wl_start : float
        lower wl
        
        wl_send : float
        upper wl"""
        
        return cls(constants.c/wl_end/1e-8, constants.c/wl_start/ 1e-8, seed=seed)
    
    def __init__(self, nu_start, nu_end, seed=250819801106):
        """Initializing photon source
        Parameters
        ----------
        
        nu_start : float
        lower nu_bin
        
        nu_end : float
        upper nu_bin
        """
        self.nu_start = nu_start
        self.nu_end = nu_end
        numpy.random.seed(seed)
        
    def __call__(self):
        # returns random frequency and random mu
        return numpy.random.uniform(self.nu_start, self.nu_end), numpy.sqrt(numpy.random.random())