from astropy import units as u

from tardis.transport.montecarlo.packet_source import (
    BlackBodySimpleSource,
)
from tardis.util.base import intensity_black_body


class BlackBodyWeightedSource(BlackBodySimpleSource):
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._reseed(
            self.base_seed
        )  # Needed because base_source doesn't seed by default

    def create_packet_nus(self, no_of_packets, l_samples=1000):
        """
        Create packet :math:`\\nu` distributed uniformly over
        bounds taken from the BlackBodySimpleSource distribution

        Parameters
        ----------
        no_of_packets : int
        l_samples : int
            number of l_samples needed for sampling from BlackBodySimpleSource

        Returns
        -------
        array of frequencies
            numpy.ndarray
        """
        nus = super().create_packet_nus(no_of_packets, l_samples)

        nu_min = nus.min()
        nu_max = nus.max()

        self.nus = (
            self.rng.uniform(nu_min.cgs.value, nu_max.cgs.value, no_of_packets)
            * nus.unit
        )

        return self.nus

    def create_packet_energies(self, no_of_packets):
        """
        Set energy weight for each packet from the relative contribution to
        the Planck Distribution

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
