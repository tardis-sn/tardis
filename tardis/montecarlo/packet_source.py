import abc

import numpy as np
import numexpr as ne
from tardis import constants as const
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)


class BasePacketSource(abc.ABC):
    def __init__(self, seed):
        self.seed = seed
        np.random.seed(seed)

    @abc.abstractmethod
    def create_packets(self, seed=None, **kwargs):
        pass

    @staticmethod
    def create_zero_limb_darkening_packet_mus(no_of_packets, rng):
        """
        Create zero-limb-darkening packet :math:`\mu` distributed
        according to :math:`\\mu=\\sqrt{z}, z \isin [0, 1]`

        Parameters
        ----------
        no_of_packets : int
            number of packets to be created
        """

        # For testing purposes
        if montecarlo_configuration.LEGACY_MODE_ENABLED:
            return np.sqrt(np.random.random(no_of_packets))
        else:
            return np.sqrt(rng.random(no_of_packets))

    @staticmethod
    def create_uniform_packet_energies(no_of_packets, rng):
        """
        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has energy of 1.

        Parameters
        ----------
        no_of_packets : int
            number of packets

        Returns
        -------
            energies for packets : numpy.ndarray
        """
        return np.ones(no_of_packets) / no_of_packets

    @staticmethod
    def create_blackbody_packet_nus(T, no_of_packets, rng, l_samples=1000):
        """
        Create packet :math:`\\nu` distributed using the algorithm described in
        Bjorkman & Wood 2001 (page 4) which references
        Carter & Cashwell 1975:
        First, generate a uniform random number, :math:`\\xi_0 \\in [0, 1]` and
        determine the minimum value of
        :math:`l, l_{\\rm min}`, that satisfies the condition

        .. math::
            \\sum_{i=1}^{l} i^{-4} \\ge {{\\pi^4}\\over{90}} m_0 \\;.

        Next obtain four additional uniform random numbers (in the range 0
        to 1) :math:`\\xi_1, \\xi_2, \\xi_3, {\\rm and } \\xi_4`.
        Finally, the packet frequency is given by

        .. math::
            x = -\\ln{(\\xi_1\\xi_2\\xi_3\\xi_4)}/l_{\\rm min}\\;.

        where :math:`x=h\\nu/kT`

        Parameters
        ----------
        T : float
            temperature
        no_of_packets : int
        l_samples : int
            number of l_samples needed in the algorithm

        Returns
        -------
            array of frequencies: numpy.ndarray
        """
        l_samples = l_samples
        l_array = np.cumsum(np.arange(1, l_samples, dtype=np.float64) ** -4)
        l_coef = np.pi**4 / 90.0

        # For testing purposes
        if montecarlo_configuration.LEGACY_MODE_ENABLED:
            xis = np.random.random((5, no_of_packets))
        else:
            xis = rng.random((5, no_of_packets))

        l = l_array.searchsorted(xis[0] * l_coef) + 1.0
        xis_prod = np.prod(xis[1:], 0)
        x = ne.evaluate("-log(xis_prod)/l")

        return x * (const.k_B.cgs.value * T) / const.h.cgs.value


class BlackBodySimpleSource(BasePacketSource):
    """
    Simple packet source that generates Blackbody packets for the Montecarlo
    part.
    """

    def create_packets(self, T, no_of_packets, rng, radius):
        """Generate black-body packet properties as arrays

        Parameters
        ----------
        T : float64
            Temperature
        no_of_packets : int
            Number of packets
        rng : numpy random number generator
        radius : float64
            Initial packet radius

        Returns
        -------
        array
            Packet radii
        array
            Packet frequencies
        array
            Packet directions
        array
            Packet energies
        """
        radii = np.ones(no_of_packets) * radius
        nus = self.create_blackbody_packet_nus(T, no_of_packets, rng)
        mus = self.create_zero_limb_darkening_packet_mus(no_of_packets, rng)
        energies = self.create_uniform_packet_energies(no_of_packets, rng)

        return radii, nus, mus, energies


class BlackBodySimpleSourceRelativistic(BlackBodySimpleSource):
    def create_packets(self, T, no_of_packets, rng, radius, time_explosion):
        """Generate relativistic black-body packet properties as arrays

        Parameters
        ----------
        T : float64
            Temperature
        no_of_packets : int
            Number of packets
        rng : numpy random number generator
        radius : float64
            Initial packet radius
        time_explosion: float64
            Time elapsed since explosion

        Returns
        -------
        array
            Packet radii
        array
            Packet frequencies
        array
            Packet directions
        array
            Packet energies
        """
        self.beta = ((radius / time_explosion) / const.c).to("")
        return super().create_packets(T, no_of_packets, rng, radius)

    def create_zero_limb_darkening_packet_mus(self, no_of_packets, rng):
        """
        Create zero-limb-darkening packet :math:`\mu^\prime` distributed
        according to :math:`\\mu^\\prime=2 \\frac{\\mu^\\prime + \\beta}{2 \\beta + 1}`.
        The complicated distribution is due to the fact that the inner boundary
        on which the packets are initialized is not comoving with the material.

        Parameters
        ----------
        no_of_packets : int
            number of packets to be created
        """
        z = rng.random(no_of_packets)
        beta = self.beta
        return -beta + np.sqrt(beta**2 + 2 * beta * z + z)

    def create_uniform_packet_energies(self, no_of_packets, rng):
        """
        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has energy of 1 multiplied by relativistic correction factors.

        Parameters
        ----------
        no_of_packets : int
            number of packets

        Returns
        -------
            energies for packets : numpy.ndarray
        """
        beta = self.beta
        gamma = 1.0 / np.sqrt(1 - beta**2)
        static_inner_boundary2cmf_factor = (2 * beta + 1) / (1 - beta**2)
        energies = np.ones(no_of_packets) / no_of_packets
        # In principle, the factor gamma should be applied to the time of
        # simulation to account for time dilation between the lab and comoving
        # frame. However, all relevant quantities (luminosities, estimators, ...)
        # are calculated as ratios of packet energies and the time of simulation.
        # Thus, we can absorb the factor gamma in the packet energies, which is
        # more convenient.
        return energies * static_inner_boundary2cmf_factor / gamma
