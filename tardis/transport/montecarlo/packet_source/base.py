import abc

import numpy as np

from tardis import constants as const
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)


class BasePacketSource(abc.ABC):
    """
    Abstract base packet source

    Parameters
    ----------
    base_seed : int
        Base Seed for random number generator
    legacy_secondary_seed : int
        Secondary seed for global numpy rng (Deprecated: Legacy reasons only)
    """

    # MAX_SEED_VAL must be multiple orders of magnitude larger than no_of_packets;
    # otherwise, each packet would not have its own seed. Here, we set the max
    # seed val to the maximum allowed by numpy.
    MAX_SEED_VAL = 2**32 - 1

    def __init__(
        self, base_seed=None, legacy_mode_enabled=False, legacy_second_seed=None
    ):
        self.base_seed = base_seed
        self.legacy_mode_enabled = legacy_mode_enabled
        if self.legacy_mode_enabled and legacy_second_seed is not None:
            np.random.seed(legacy_second_seed)
        else:
            np.random.seed(self.base_seed)

    def _reseed(self, seed):
        self.rng = np.random.default_rng(seed=seed)

    # One should either implement create_packet_velocities or create_packet_radii
    def create_packet_radii(self, no_of_packets, *args, **kwargs):
        pass

    def create_packet_velocities(self, no_of_packets, *args, **kwargs):
        pass

    @abc.abstractmethod
    def create_packet_nus(self, no_of_packets, *args, **kwargs):
        pass

    @abc.abstractmethod
    def create_packet_mus(self, no_of_packets, *args, **kwargs):
        pass

    @abc.abstractmethod
    def create_packet_energies(self, no_of_packets, *args, **kwargs):
        pass

    def create_packets(self, no_of_packets, seed_offset=0, *args, **kwargs):
        """Generate packet properties as arrays

        Parameters
        ----------
        no_of_packets : int
            Number of packets

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
        # the iteration (passed as seed_offset) is added each time to preserve randomness
        # across different simulations with the same temperature,
        # for example. We seed the random module instead of the numpy module
        # because we call random.sample, which references a different internal
        # state than in the numpy.random module.
        self._reseed(self.base_seed + seed_offset)
        packet_seeds = self.rng.choice(
            self.MAX_SEED_VAL, no_of_packets, replace=True
        )

        radii = self.create_packet_radii(no_of_packets, *args, **kwargs).value
        nus = self.create_packet_nus(no_of_packets, *args, **kwargs).value
        mus = self.create_packet_mus(no_of_packets, *args, **kwargs)
        energies = self.create_packet_energies(
            no_of_packets, *args, **kwargs
        ).value
        # Check if all arrays have the same length
        assert (
            len(radii) == len(nus) == len(mus) == len(energies) == no_of_packets
        )
        radiation_field_luminosity = self.calculate_radfield_luminosity().value
        return PacketCollection(
            radii,
            nus,
            mus,
            energies,
            packet_seeds,
            radiation_field_luminosity,
        )

    def calculate_radfield_luminosity(self):
        """
        Calculate inner luminosity.

        Parameters
        ----------
        model : model.SimulationState

        Returns
        -------
        astropy.units.Quantity
        """
        return (
            4
            * np.pi
            * const.sigma_sb
            * self.radius**2
            * self.temperature**4
        ).to("erg/s")




