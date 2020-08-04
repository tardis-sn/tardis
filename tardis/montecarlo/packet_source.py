import abc

import numpy as np
import numexpr as ne
from tardis import constants as const
from tardis.montecarlo import montecarlo_configuration as mc_config_module
from numba import njit, float64



class BasePacketSource(abc.ABC):

    def __init__(self, seed):
        self.seed = seed
        np.random.seed(seed)
        
    @abc.abstractmethod
    def create_packets(self, seed=None, **kwargs):
        pass

    @staticmethod
    @njit
    def create_zero_limb_darkening_packet_mus(seed):
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

    @staticmethod
    def create_uniform_packet_energies(no_of_packets):
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


    @staticmethod
    def create_blackbody_packet_nus(T, xis, l_samples=1000):
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
        l_array = np.cumsum(np.arange(1, l_samples, dtype=np.float64)**-4)
        l_coef = np.pi**4 / 90.0

        l = np.searchsorted(l_array, xis[0]*l_coef) + 1.
        xis_prod = np.prod(xis[1:], 0)
        x = ne.evaluate('-log(xis_prod)/l')

        return x * (const.k_B.cgs.value * T) / const.h.cgs.value


class BlackBodySimpleSource(BasePacketSource):
    """
    Simple packet source that generates Blackbody packets for the Montecarlo 
    part.
    """

    def create_packets(self, T, no_of_packets):
        self.T = T
        self.seeds = mc_config_module.packet_seeds
        nus = self.create_nus()
        mus = self.create_mus()
        energies = self.create_uniform_packet_energies(no_of_packets)
        return nus, mus, energies

    def create_nus(self):
        xis = np.zeros((5, len(self.seeds)))
        for i, seed in enumerate(self.seeds):
            xis[:, i] = self.random_seeded_array(seed)
        return self.create_blackbody_packet_nus(self.T, xis)

    @staticmethod
    @njit
    def random_seeded_array(seed):
        np.random.seed(seed)
        return np.random.random(5)

    def create_mus(self):
        return np.array([self.create_zero_limb_darkening_packet_mus(seed)
                         for seed in self.seeds])
