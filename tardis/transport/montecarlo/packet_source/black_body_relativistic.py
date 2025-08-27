from tardis import constants as const
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.transport.montecarlo.packet_source.black_body import BlackBodySimpleSource


import numpy as np
from astropy import units as u


class BlackBodySimpleSourceRelativistic(BlackBodySimpleSource, HDFWriterMixin):
    """
    Simple packet source that generates Blackbody packets for the Montecarlo
    part.

    Parameters
    ----------
    time_explosion : astropy.units.Quantity
        Time elapsed since explosion
    radius : astropy.units.Quantity
        Initial packet radius
    temperature : astropy.units.Quantity
        Absolute Temperature.
    base_seed : int
        Base Seed for random number generator
    legacy_secondary_seed : int
        Secondary seed for global numpy rng (Deprecated: Legacy reasons only)
    """

    hdf_properties = ["time_explosion", "radius", "temperature", "base_seed"]

    @classmethod
    def from_simulation_state(cls, simulation_state, *args, **kwargs):
        return cls(
            simulation_state.time_explosion,
            simulation_state.r_inner[0],
            simulation_state.t_inner,
            *args,
            **kwargs,
        )

    def __init__(self, time_explosion=None, **kwargs):
        self.time_explosion = time_explosion
        super().__init__(**kwargs)

    def create_packets(self, no_of_packets, *args, **kwargs):
        """Generate relativistic black-body packet properties as arrays

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
        if self.radius is None or self.time_explosion is None:
            raise ValueError("Black body Radius or Time of Explosion isn't set")
        self.beta = ((self.radius / self.time_explosion) / const.c).to(u.dimensionless_unscaled).value
        return super().create_packets(no_of_packets, *args, **kwargs)

    def create_packet_mus(self, no_of_packets):
        """
        Create zero-limb-darkening packet :math:`\\mu^\\prime` distributed
        according to :math:`\\mu^\\prime=2 \\frac{\\mu^\\prime + \\beta}{2 \\beta + 1}`.
        The complicated distribution is due to the fact that the inner boundary
        on which the packets are initialized is not comoving with the material.

        Parameters
        ----------
        no_of_packets : int
            number of packets to be created

        Returns
        -------
        Directions for packets
            numpy.ndarray
        """
        z = self.rng.random(no_of_packets)
        beta = self.beta
        return -beta + np.sqrt(beta**2 + 2 * beta * z + z)

    def create_packet_energies(self, no_of_packets):
        """
        Uniformly distribute energy in arbitrary units where the ensemble of
        packets has energy of 1 multiplied by relativistic correction factors.

        Parameters
        ----------
        no_of_packets : int
            number of packets

        Returns
        -------
        energies for packets
            numpy.ndarray
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