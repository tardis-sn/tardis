#creating photons

import numpy as np

from astropy import units

import logging
from util import intensity_black_body

logger = logging.getLogger(__name__)

class SimplePacketSource:
    """Initializing photon source
        Parameters
        ----------

        nu_start : float
            lowest frequency

        nu_end : float
            highest_frequency
    """

    @classmethod
    def from_wavelength(cls, wavelength_start, wavelength_end,  seed=250819801106, blackbody_sampling=int(1e6)):
        """Initializing from wavelength

        Parameters
        ----------
        
        wavelength_start : `float`
            start of the wavelength
        
        wavelength_end : `float`
        upper wl"""

        nu_start = wavelength_end.to('Hz', units.spectral()).value
        nu_end = wavelength_start.to('Hz', units.spectral()).value

        return cls(nu_start, nu_end, seed=seed, blackbody_sampling=blackbody_sampling)

    def __init__(self, nu_start, nu_end, seed=250819801106, blackbody_sampling=int(1e6)):
        self.nu_start = nu_start
        self.nu_end = nu_end
        self.blackbody_sampling = blackbody_sampling
        np.random.seed(seed)


    def create_packets(self, number_of_packets, t_rad, seed=None):
        """
        Creating a new random number of packets, with a certain temperature

        Parameters
        ----------

        number_of_packets : any number
            number of packets

        t_rad : `float`
            radiation temperature

        """
        if seed is not None:
            np.random.seed(seed)

        number_of_packets = int(number_of_packets)


        self.packet_nus = self.random_blackbody_nu(t_rad, number_of_packets)

        self.packet_mus = np.sqrt(np.random.random(size=number_of_packets))
        self.packet_energies = np.ones(number_of_packets) / number_of_packets


    def random_blackbody_nu(self, T, number_of_packets):
        """
        Creating the random nus for the energy packets

        Parameters
        ----------

        T : `float`
            temperature of the blackbody

        number_of_packets : `int`
            the number of packets
        """
        logger.info('Calculating %d packets for t_inner=%.2f', number_of_packets, T)
        nu = np.linspace(self.nu_start, self.nu_end, num=self.blackbody_sampling)
        intensity = intensity_black_body(nu, T)
        cum_blackbody = np.cumsum(intensity)
        norm_cum_blackbody = cum_blackbody / cum_blackbody.max()
        return nu[norm_cum_blackbody.searchsorted(np.random.random(number_of_packets))] + \
               np.random.random(size=number_of_packets) * (nu[1] - nu[0])


