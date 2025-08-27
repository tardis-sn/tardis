

from typing import Any

import numexpr as ne
import numpy as np
from astropy import units as u

from tardis import constants as const
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.transport.montecarlo.packet_source.base import BasePacketSource
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)


class BlackBodySimpleSource(BasePacketSource, HDFWriterMixin):
    """
    Simple packet source that generates blackbody packets for Monte Carlo simulations.

    This class creates packets with properties derived from blackbody radiation,
    including appropriate frequency distribution, uniform radii, and cosine-weighted
    direction distribution.

    Parameters
    ----------
    radius : astropy.units.Quantity, optional
        Initial packet radius. Default is None.
    temperature : astropy.units.Quantity, optional
        Blackbody temperature. Default is None.
    **kwargs
        Additional keyword arguments passed to the parent class.

    Attributes
    ----------
    radius : astropy.units.Quantity
        Initial packet radius.
    temperature : astropy.units.Quantity
        Blackbody temperature.
    """

    hdf_properties = ["radius", "temperature", "base_seed"]
    hdf_name = "black_body_simple_source"

    @classmethod
    def from_simulation_state(
        cls, simulation_state, *args: Any, **kwargs: Any
    ) -> "BlackBodySimpleSource":
        """
        Create BlackBodySimpleSource from simulation state.

        Parameters
        ----------
        simulation_state : SimulationState
            The simulation state object containing inner radius and temperature.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        BlackBodySimpleSource
            New instance initialized with simulation state parameters.
        """
        return cls(
            simulation_state.r_inner[0],
            simulation_state.t_inner,
            *args,
            **kwargs,
        )

    def __init__(
        self,
        radius: "u.Quantity | None" = None,
        temperature: "u.Quantity | None" = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize BlackBodySimpleSource.

        Parameters
        ----------
        radius : astropy.units.Quantity, optional
            Initial packet radius. Default is None.
        temperature : astropy.units.Quantity, optional
            Absolute temperature. Default is None.
        **kwargs : Any
            Additional keyword arguments passed to parent class.
        """
        self.radius = radius
        self.temperature = temperature
        super().__init__(**kwargs)

    def create_packets(self, no_of_packets: int, *args: Any, **kwargs: Any) -> PacketCollection:
        """
        Create packet collection.

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
            Collection of packets.

        Raises
        ------
        ValueError
            If radius or temperature is not set.
        """
        if self.radius is None or self.temperature is None:
            raise ValueError("Black body Radius or Temperature isn't set")
        return super().create_packets(no_of_packets, *args, **kwargs)

    def create_packet_radii(self, no_of_packets: int) -> u.Quantity:
        """
        Create packet radii.

        All packets are created at the same radius (inner boundary).

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.

        Returns
        -------
        astropy.units.Quantity
            Array of packet radii in CGS units.
        """
        return np.ones(no_of_packets) * self.radius.cgs

    def create_packet_nus(self, no_of_packets: int, l_samples: int = 1000) -> u.Quantity:
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
        no_of_packets : int
        l_samples : int
            number of l_samples needed in the algorithm

        Returns
        -------
        array of frequencies
            numpy.ndarray
        """
        l_array = np.cumsum(np.arange(1, l_samples, dtype=np.float64) ** -4)
        l_coef = np.pi**4 / 90.0

        # For testing purposes
        if self.legacy_mode_enabled:
            xis = np.random.random((5, no_of_packets))
        else:
            xis = self.rng.random((5, no_of_packets))

        l = l_array.searchsorted(xis[0] * l_coef) + 1.0
        xis_prod = np.prod(xis[1:], 0)
        x = ne.evaluate("-log(xis_prod)/l")

        return (x * (const.k_B * self.temperature) / const.h).cgs

    def create_packet_mus(self, no_of_packets: int) -> np.ndarray:
        """
        Create zero-limb-darkening packet direction cosines.

        Direction cosines are distributed according to :math:`\\mu=\\sqrt{z}`,
        where :math:`z \\in [0, 1]` is uniformly distributed.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.

        Returns
        -------
        numpy.ndarray
            Array of direction cosines for packets.
        """
        # For testing purposes
        if self.legacy_mode_enabled:
            return np.sqrt(np.random.random(no_of_packets))
        return np.sqrt(self.rng.random(no_of_packets))

    def create_packet_energies(self, no_of_packets: int) -> u.Quantity:
        """
        Create packet energies with uniform distribution.

        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has total energy of 1 erg.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create.

        Returns
        -------
        astropy.units.Quantity
            Array of packet energies in erg.
        """
        return np.ones(no_of_packets) / no_of_packets * u.erg

    def set_temperature_from_luminosity(self, luminosity: u.Quantity) -> None:
        """
        Set blackbody packet source temperature from luminosity.

        Uses the Stefan-Boltzmann law to derive temperature from the given
        luminosity and the source radius.

        Parameters
        ----------
        luminosity : astropy.units.Quantity
            Total luminosity to match.
        """
        self.temperature = (
            (luminosity / (4 * np.pi * self.radius**2 * const.sigma_sb))
            ** 0.25
        ).to("K")
