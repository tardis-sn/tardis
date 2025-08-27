from typing import Any

import numpy as np
from astropy import units as u

from tardis import constants as const
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.transport.montecarlo.packet_source.black_body import (
    BlackBodySimpleSource,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)


class BlackBodySimpleSourceRelativistic(BlackBodySimpleSource, HDFWriterMixin):
    """
    Relativistic blackbody packet source for Monte Carlo simulations.

    This class generates blackbody packets with relativistic corrections
    for sources that are moving with respect to the lab frame. It accounts
    for the motion of the inner boundary where packets are created.

    Parameters
    ----------
    time_explosion : astropy.units.Quantity
        Time elapsed since explosion.
    radius : astropy.units.Quantity
        Initial packet radius.
    temperature : astropy.units.Quantity
        Absolute temperature.
    base_seed : int, optional
        Base seed for random number generator.
    legacy_secondary_seed : int, optional
        Secondary seed for global numpy rng (Deprecated: Legacy reasons only).

    Attributes
    ----------
    time_explosion : astropy.units.Quantity
        Time elapsed since explosion.
    beta : float
        Velocity of the inner boundary as a fraction of speed of light.
    """

    hdf_properties = ["time_explosion", "radius", "temperature", "base_seed"]

    @classmethod
    def from_simulation_state(
        cls, simulation_state, *args: Any, **kwargs: Any
    ) -> "BlackBodySimpleSourceRelativistic":
        """
        Create BlackBodySimpleSourceRelativistic from simulation state.

        Parameters
        ----------
        simulation_state : SimulationState
            The simulation state object containing explosion time, inner radius, and temperature.
        *args : Any
            Additional positional arguments.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        BlackBodySimpleSourceRelativistic
            New instance initialized with simulation state parameters.
        """
        return cls(
            simulation_state.time_explosion,
            simulation_state.r_inner[0],
            simulation_state.t_inner,
            *args,
            **kwargs,
        )

    def __init__(
        self, time_explosion: "u.Quantity | None" = None, **kwargs: Any
    ) -> None:
        """
        Initialize BlackBodySimpleSourceRelativistic.

        Parameters
        ----------
        time_explosion : astropy.units.Quantity, optional
            Time elapsed since explosion. Default is None.
        **kwargs : Any
            Additional keyword arguments passed to parent class.
        """
        self.time_explosion = time_explosion
        super().__init__(**kwargs)

    def create_packets(
        self, no_of_packets: int, *args: Any, **kwargs: Any
    ) -> PacketCollection:
        """
        Generate relativistic black-body packet properties as arrays.

        Calculates the velocity (beta) of the inner boundary and applies
        relativistic corrections to packet creation.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.
        *args : Any
            Additional positional arguments.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        PacketCollection
            Collection of packets with relativistic corrections applied.

        Raises
        ------
        ValueError
            If radius or time_explosion is not set.
        """
        if self.radius is None or self.time_explosion is None:
            raise ValueError("Black body Radius or Time of Explosion isn't set")
        self.beta = ((self.radius / self.time_explosion) / const.c).to(u.dimensionless_unscaled).value
        return super().create_packets(no_of_packets, *args, **kwargs)

    def create_packet_mus(self, no_of_packets: int) -> np.ndarray:
        """
        Create relativistic packet direction cosines.

        Direction cosines are distributed according to the relativistic transformation
        :math:`\\mu^\\prime=2 \\frac{\\mu^\\prime + \\beta}{2 \\beta + 1}`.
        The distribution accounts for the fact that the inner boundary
        on which the packets are initialized is not comoving with the material.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.

        Returns
        -------
        numpy.ndarray
            Array of relativistically corrected direction cosines for packets.
        """
        z = self.rng.random(no_of_packets)
        beta = self.beta
        return -beta + np.sqrt(beta**2 + 2 * beta * z + z)

    def create_packet_energies(self, no_of_packets: int) -> u.Quantity:
        """
        Create relativistic packet energies with uniform distribution.

        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has total energy corrected for relativistic effects. Applies
        corrections for the static inner boundary to comoving frame transformation
        and time dilation.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.

        Returns
        -------
        astropy.units.Quantity
            Array of packet energies in erg with relativistic corrections.
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
        return energies * static_inner_boundary2cmf_factor / gamma * u.erg
