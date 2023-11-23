import numpy as np

from astropy import constants as const, units as u
from scipy.special import zeta

DILUTION_FACTOR_ESTIMATOR_CONSTANT = (
    (const.c**2 / (2 * const.h))
    * (15 / np.pi**4)
    * (const.h / const.k_B) ** 4
    / (4 * np.pi)
).cgs.value

T_RADIATIVE_ESTIMATOR_CONSTANT = (
    (np.pi**4 / (15 * 24 * zeta(5, 1))) * (const.h / const.k_B)
).cgs.value


class MonteCarloTransportState:
    def __init__(self, packet_collection, estimators, volume):
        self.packet_collection = packet_collection
        self.estimators = estimators
        self.volume = volume

    def calculate_radiationfield_properties(self):
        """
        Calculate an updated radiation field from the :math:
        `\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation.
        The details of the calculation can be found in the documentation.

        Parameters
        ----------
        nubar_estimator : np.ndarray (float)
        j_estimator : np.ndarray (float)

        Returns
        -------
        t_radiative : astropy.units.Quantity (float)
        dilution_factor : numpy.ndarray (float)
        """

        estimated_t_radiative = (
            T_RADIATIVE_ESTIMATOR_CONSTANT
            * self.estimators.nu_bar_estimator
            / self.estimators.j_estimator
        ) * u.K
        dilution_factor = self.estimators.j_estimator / (
            4
            * const.sigma_sb.cgs.value
            * estimated_t_radiative.value**4
            * (self.packet_collection.time_of_simulation)
            * self.volume.cgs.value
        )

        return estimated_t_radiative, dilution_factor

    @property
    def packet_luminosity(self):
        return (
            (
                self.packet_collection.output_energies
                / self.packet_collection.time_of_simulation
            )
            * u.erg
            / u.s
        )

    @property
    def emitted_packet_mask(self):
        return self.packet_collection.output_energies >= 0

    @property
    def emitted_packet_nu(self):
        return (
            self.packet_collection.output_nus[self.emitted_packet_mask] * u.Hz
        )

    @property
    def reabsorbed_packet_nu(self):
        return (
            self.packet_collection.output_nus[~self.emitted_packet_mask] * u.Hz
        )

    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]

    def calculate_emitted_luminosity(
        self, luminosity_nu_start, luminosity_nu_end
    ):
        """
        Calculate emitted luminosity.

        Parameters
        ----------
        luminosity_nu_start : astropy.units.Quantity
        luminosity_nu_end : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        luminosity_wavelength_filter = (
            self.emitted_packet_nu > luminosity_nu_start
        ) & (self.emitted_packet_nu < luminosity_nu_end)

        return self.emitted_packet_luminosity[
            luminosity_wavelength_filter
        ].sum()

    def calculate_reabsorbed_luminosity(
        self, luminosity_nu_start, luminosity_nu_end
    ):
        """
        Calculate reabsorbed luminosity.

        Parameters
        ----------
        luminosity_nu_start : astropy.units.Quantity
        luminosity_nu_end : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        luminosity_wavelength_filter = (
            self.reabsorbed_packet_nu > luminosity_nu_start
        ) & (self.reabsorbed_packet_nu < luminosity_nu_end)

        return self.reabsorbed_packet_luminosity[
            luminosity_wavelength_filter
        ].sum()
