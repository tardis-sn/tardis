import numexpr as ne
import numpy as np
from astropy import constants as const
from astropy import units as u

from tardis.transport.montecarlo.packet_source import (
    BasePacketSource,
    HDFWriterMixin,
)
from tardis.util.base import intensity_black_body


class BlackBodyWeightedSource(BasePacketSource, HDFWriterMixin):
    """
    Simple packet source that generates Blackbody packets for the Montecarlo
    part.

    Parameters
    ----------
    radius : astropy.units.Quantity
        Initial packet radius
    temperature : astropy.units.Quantity
        Absolute Temperature.
    base_seed : int
        Base Seed for random number generator
    legacy_secondary_seed : int
        Secondary seed for global numpy rng (Deprecated: Legacy reasons only)
    """

    hdf_properties = ["radius", "temperature", "base_seed"]
    hdf_name = "black_body_weighted_source"

    @classmethod
    def from_simulation_state(cls, simulation_state, *args, **kwargs):
        return cls(
            simulation_state.r_inner[0],
            simulation_state.t_inner,
            *args,
            **kwargs,
        )

    def __init__(self, radius=None, temperature=None, **kwargs):
        self.radius = radius
        self.temperature = temperature
        super().__init__(**kwargs)

    def create_packets(self, no_of_packets, *args, **kwargs):
        if self.radius is None or self.temperature is None:
            raise ValueError("Black body Radius or Temperature isn't set")
        return super().create_packets(no_of_packets, *args, **kwargs)

    def create_packet_radii(self, no_of_packets):
        """
        Create packet radii

        Parameters
        ----------
        no_of_packets : int
            number of packets to be created

        Returns
        -------
        Radii for packets
            numpy.ndarray
        """
        return np.ones(no_of_packets) * self.radius.cgs

    def create_packet_nus(self, no_of_packets, l_samples=1000):
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

        nus = (x * (const.k_B * self.temperature) / const.h).cgs

        nu_min = nus.min()
        nu_max = nus.max()

        self.nus = np.random.uniform(nu_min.cgs.value, nu_max.cgs.value, no_of_packets)*nus.unit

        return self.nus

    def create_packet_mus(self, no_of_packets):
        """
        Create zero-limb-darkening packet :math:`\\mu` distributed
        according to :math:`\\mu=\\sqrt{z}, z \\isin [0, 1]`

        Parameters
        ----------
        no_of_packets : int
            number of packets to be created

        Returns
        -------
        Directions for packets
            numpy.ndarray
        """
        # For testing purposes
        if self.legacy_mode_enabled:
            return np.sqrt(np.random.random(no_of_packets))
        else:
            return np.sqrt(self.rng.random(no_of_packets))

    def create_packet_energies(self, no_of_packets):
        """
        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has energy of 1.

        Parameters
        ----------
        no_of_packets : int
            number of packets

        Returns
        -------
        energies for packets
            numpy.ndarray
        """
        try:
            self.nus
        except AttributeError:
            self.nus = self.create_packet_nus(no_of_packets)

        intensity = intensity_black_body(self.nus.cgs.value, self.temperature)
        return intensity / intensity.sum() * u.erg

    def set_temperature_from_luminosity(self, luminosity: u.Quantity):
        """
        Set blackbody packet source temperature from luminosity

        Parameters
        ----------
        luminosity : u.Quantity

        """
        self.temperature = (
            (luminosity / (4 * np.pi * self.radius**2 * const.sigma_sb))
            ** 0.25
        ).to("K")

