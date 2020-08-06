import abc

import numpy as np
import numexpr as ne
from tardis import constants as const
from numba import jitclass, njit, int64, float64

BOLTZMANN_CONSTANT = const.k_B.cgs.value

PLANCK_CONSTANT = const.h.cgs.value

packet_source_spec = [
    ('seed', int64),
]

# @jitclass(packet_source_spec)
class BasePacketSource(object):

    def __init__(self, seed):
        self.seed = seed
        np.random.seed(seed)
        self
        
    # @abc.abstractmethod
    # def create_packets(self, **kwargs):
    #     pass

    def create_zero_limb_darkening_packet_mus(self, seed):
        """
        Create a zero-limb-darkening packet :math:`\mu` distributed
        according to :math:`\\mu=\\sqrt{z}, z \isin [0, 1]`
        
        Parameters
        ----------
        seed : int
            value to seed random number generator.
        """
        np.random.seed(seed)
        return np.sqrt(np.random.random())

    def create_uniform_packet_energies(self, no_of_packets):
        """
        Uniformly distribute energy in arbitrary units where the ensemble of 
        packets has energy of 1. 

        
        Parameters
        ----------
        no_of_packets : int
            number of packets
        
        Returns
        -------
            : numpy.ndarray
            energies for packets
        """
        return np.ones(no_of_packets) / no_of_packets

    def create_blackbody_packet_nus(self, T, xis, l_samples=1000.):
        """

        Create packet :math:`\\nu` distributed using the algorithm described in 
        Bjorkman & Wood 2001 (page 4) which references 
        Carter & Cashwell 1975:

        First, generate a uniform random number, :math:`\\xi_0 \\in [0, 1]` and
        determine the minimum value of 
        :math:`l, l_{\\rm min}`, that satisfies the condition

        .. math::
            \\sum_{i=1}^{l} i^{-4} \\ge {{\\pi^4}\\over{90}} m_0 \\; .
        
        Next obtain four additional uniform random numbers (in the range 0 
        to 1) :math:`\\xi_1, \\xi_2, \\xi_3, {\\rm and } \\xi_4`.  
        
        Finally, the packet frequency is given by
        
        .. math::
            x = -\\ln{(\\xi_1\\xi_2\\xi_3\\xi_4)}/l_{\\rm min}\\; .


        where :math:`x=h\\nu/kT`

        Parameters
        ----------
        T : float
            temperature
        seed: int
        l_samples: int
            number of l_samples needed in the algorithm

        Returns
        -------

            : numpy.ndarray
            array of frequencies
        """
        l_samples = l_samples
        l_array = np.cumsum(np.arange(1., l_samples, dtype=float64)**-4)
        l_coef = np.pi**4 / 90.0

        l = np.searchsorted(l_array, xis[0]*l_coef) + 1.
        xis_prod = np.prod(xis[1:])
        x = -np.log(xis_prod)/l

        return x * (BOLTZMANN_CONSTANT * T) / PLANCK_CONSTANT

@jitclass(packet_source_spec)
class BlackBodySimpleSource(BasePacketSource):
    """
    Simple packet source that generates Blackbody packets for the Montecarlo 
    part.
    """

    def create_packets(self, T, no_of_packets, seeds):
        """
        Creates a number of r_packet properties as sampled from a blackbody.

        Inputs:
            :T: (float) temperature of blackbody to be modeled.
            :no_of_packets: (int) number of Monte Carlo packets to model.

        Outputs:
            :nus: (1D array) array of frequencies associated with packets.
        """
        # loop through: set seed, create a mu, create a nu
        nus, mus, energies = np.empty(no_of_packets), np.empty(no_of_packets), \
                                np.empty(no_of_packets)
        for i, seed in enumerate(seeds):
            nus[i], mus[i], energies[i] = self.random_packet_properties(seed,
                                                                        T,
                                                                        self.create_blackbody_packet_nus,
                                                                        no_of_packets)
        return nus, mus, energies

    def random_packet_properties(self, seed, T, blackbody_sampler, no_of_packets):
        """
        Created random packet properties (energy, frequency, mu) given a seed.
        The blackbody_sampler is passed directly (with this method being a
        staticmethod) to circumvent having to @jitclass BasePacketSource or
        BlackBodySimpleSource.

        Inputs:
            :seed: (int) value to seed the random number generator. Should be
                    distinct for each packet.
            :T: (float) temperature at which the blackbody should be modeled.
            :blackbody_sampler: (function) function to sample the blackbody.
            :no_of_packets: (int) total number of r_packets being simulated.

        Outputs:
            :nu: frequency of r_packet with this seed.
            :mu: mu of r_packet with this seed.
            :energy: energy of r_packet with this seed.
        """

        np.random.seed(seed)

        xis = np.random.random(5)
        nu = blackbody_sampler(T, xis)

        mu = np.sqrt(np.random.random())  # zero limb darkening

        energy = 1/no_of_packets  # uniform packet energies
        return nu, mu, energy

