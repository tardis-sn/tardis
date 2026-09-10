import abc
from typing import Any

import numpy as np
from astropy import units as u

from tardis import constants as const
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)


class BasePacketSource(abc.ABC):
    """
    Abstract base packet source.

    This abstract base class defines the interface for packet sources used in
    TARDIS Monte Carlo radiative transfer. Packet sources are responsible
    for creating radiation packets with specific properties.

    Parameters
    ----------
    base_seed : int, optional
        Base seed for random number generator. Default is None.
    legacy_mode_enabled : bool, optional
        Whether to enable legacy mode for compatibility. Default is False.
    legacy_second_seed : int, optional
        Secondary seed for global numpy rng (Deprecated: Legacy reasons only).
        Default is None.

    Attributes
    ----------
    MAX_SEED_VAL : int
        Maximum seed value allowed by numpy (2**32 - 1).
    base_seed : int or None
        Base seed for random number generator.
    legacy_mode_enabled : bool
        Whether legacy mode is enabled.
    rng : numpy.random.Generator
        Random number generator instance.
    """

    # MAX_SEED_VAL must be multiple orders of magnitude larger than no_of_packets;
    # otherwise, each packet would not have its own seed. Here, we set the max
    # seed val to the maximum allowed by numpy.
    MAX_SEED_VAL: int = 2**32 - 1

    def __init__(
        self,
        base_seed: "int | None" = None,
        legacy_mode_enabled: bool = False,
        legacy_second_seed: "int | None" = None,
    ) -> None:
        self.base_seed = base_seed
        self.legacy_mode_enabled = legacy_mode_enabled
        if self.legacy_mode_enabled and legacy_second_seed is not None:
            np.random.seed(legacy_second_seed)
        else:
            np.random.seed(self.base_seed)

    def _reseed(self, seed: int) -> None:
        """
        Reseed the random number generator.

        Parameters
        ----------
        seed : int
            Seed value for the random number generator.
        """
        self.rng = np.random.default_rng(seed=seed)

    # One should either implement create_packet_velocities or create_packet_radii
    def create_packet_radii(self, no_of_packets: int, *args: Any, **kwargs: Any):
        """
        Create packet radii.

        This method should be implemented by subclasses that create packets
        with specific radii. Either this method or create_packet_velocities
        should be implemented.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        array-like
            Packet radii values.

        Raises
        ------
        NotImplementedError
            If the method is not implemented by the subclass.
        """
        raise NotImplementedError(
            "Subclasses must implement either create_packet_radii or create_packet_velocities"
        )

    def create_packet_velocities(self, no_of_packets: int, *args: Any, **kwargs: Any):
        """
        Create packet velocities.

        This method should be implemented by subclasses that create packets
        with specific velocities. Either this method or create_packet_radii
        should be implemented.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        array-like
            Packet velocity values.

        Raises
        ------
        NotImplementedError
            If the method is not implemented by the subclass.
        """
        raise NotImplementedError(
            "Subclasses must implement either create_packet_radii or create_packet_velocities"
        )

    @abc.abstractmethod
    def create_packet_nus(self, no_of_packets: int, *args: Any, **kwargs: Any):
        """
        Create packet frequencies.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        array-like
            Packet frequency values.
        """

    @abc.abstractmethod
    def create_packet_mus(self, no_of_packets: int, *args: Any, **kwargs: Any):
        """
        Create packet direction cosines.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        array-like
            Packet direction cosine values.
        """

    @abc.abstractmethod
    def create_packet_energies(self, no_of_packets: int, *args: Any, **kwargs: Any):
        """
        Create packet energies.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        array-like
            Packet energy values.
        """

    def create_packets(
        self,
        no_of_packets: int,
        seed_offset: int = 0,
        *args: Any,
        **kwargs: Any,
    ) -> PacketCollection:
        """
        Generate packet properties as arrays.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        seed_offset : int, optional
            Offset added to the base seed for randomness across iterations.
            Default is 0.
        *args
            Additional positional arguments passed to packet creation methods.
        **kwargs
            Additional keyword arguments passed to packet creation methods.

        Returns
        -------
        PacketCollection
            Collection containing packet radii, frequencies, directions,
            energies, seeds, and radiation field luminosity.
        """
        # the iteration (passed as seed_offset) is added each time to preserve randomness
        # across different simulations with the same temperature,
        # for example. We seed the random module instead of the numpy module
        # because we call random.sample, which references a different internal
        # state than in the numpy.random module.
        if self.base_seed is None:
            raise ValueError("base_seed must be set before creating packets")
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

    def calculate_radfield_luminosity(self) -> u.Quantity:
        """
        Calculate inner luminosity from blackbody radiation.

        Uses the Stefan-Boltzmann law to calculate the luminosity from
        the inner boundary radius and temperature.

        Returns
        -------
        astropy.units.Quantity
            Inner luminosity in erg/s.
        """
        return (
            4
            * np.pi
            * const.sigma_sb
            * self.radius**2
            * self.temperature**4
        ).to("erg/s")




