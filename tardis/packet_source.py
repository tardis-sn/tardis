#creating photons

import numpy as np

from astropy import units

import plasma



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
    def from_wavelength(cls, wavelength_start, wavelength_end, wavelength_unit='angstrom', seed=250819801106,
                        blackbody_sampling=int(1e6)):


        """Initializing from wavelength

        Parameters
        ----------
        
        wavelength_start : `float`
            start of the wavelength
        
        wavelength_end : `float`
        upper wl"""

        nu_start = units.Unit(wavelength_unit).to('Hz', wavelength_end, units.spectral())
        nu_end = units.Unit(wavelength_unit).to('Hz', wavelength_start, units.spectral())

        return cls(nu_start, nu_end, seed=seed)

    def __init__(self, nu_start, nu_end, seed=250819801106, blackbody_sampling=int(1e6)):

        self.nu_start = nu_start
        self.nu_end = nu_end
        self.blackbody_sampling = blackbody_sampling
        np.random.seed(seed)


    def create_packets(self, number_of_packets, T):
        """
        Creating a new random number of packets, with a certain temperature

        Parameters
        ----------

        number_of_packets : any number
            number of packets

        T : `float`
            temperature

        """
        number_of_packets = int(number_of_packets)
        self.packet_nus = self.random_blackbody_nu(T, number_of_packets)
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
        nu = np.linspace(self.nu_start, self.nu_end, num=self.blackbody_sampling)
        intensity = plasma.intensity_black_body(nu, T)
        cum_blackbody = np.cumsum(intensity)
        norm_cum_blackbody = cum_blackbody / cum_blackbody.max()
        return nu[norm_cum_blackbody.searchsorted(np.random.random(number_of_packets))] + \
               np.random.random(size=number_of_packets) * (nu[1] - nu[0])


