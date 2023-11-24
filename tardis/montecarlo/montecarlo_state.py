import numpy as np
import warnings

from astropy import units as u
from tardis import constants as const
from scipy.special import zeta

from tardis.montecarlo.spectrum import TARDISSpectrum

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
    def __init__(
        self, packet_collection, estimators, volume, spectrum_frequency
    ):
        self.packet_collection = packet_collection
        self.estimators = estimators
        self.volume = volume
        self.spectrum_frequency = spectrum_frequency

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

    @property
    def montecarlo_reabsorbed_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.reabsorbed_packet_nu,
                weights=self.reabsorbed_packet_luminosity,
                bins=self.spectrum_frequency,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_emitted_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.emitted_packet_nu,
                weights=self.emitted_packet_luminosity,
                bins=self.spectrum_frequency,
            )[0],
            "erg / s",
        )

    @property
    def spectrum(self):
        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_emitted_luminosity
        )

    @property
    def spectrum_reabsorbed(self):
        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_reabsorbed_luminosity
        )

    @property
    def spectrum_virtual(self):
        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                "MontecarloTransport.spectrum_virtual"
                "is zero. Please run the montecarlo simulation with"
                "no_of_virtual_packets > 0",
                UserWarning,
            )

        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_virtual_luminosity
        )

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None:
            # This was changed from unpacking to specific attributes as compute
            # is not used in calculate_spectrum
            self._spectrum_integrated = self.integrator.calculate_spectrum(
                self.spectrum_frequency[:-1],
                points=self.integrator_settings.points,
                interpolate_shells=self.integrator_settings.interpolate_shells,
            )
        return self._spectrum_integrated

    @property
    def integrator(self):
        if self._integrator is None:
            warnings.warn(
                "MontecarloTransport.integrator: "
                "The FormalIntegrator is not yet available."
                "Please run the montecarlo simulation at least once.",
                UserWarning,
            )
        if self.enable_full_relativity:
            raise NotImplementedError(
                "The FormalIntegrator is not yet implemented for the full "
                "relativity mode. "
                "Please run with config option enable_full_relativity: "
                "False."
            )
        return self._integrator

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

    @property
    def output_nu(self):
        return u.Quantity(self._output_nu, u.Hz)

    @property
    def output_energy(self):
        return u.Quantity(self._output_energy, u.erg)

    @property
    def virtual_packet_nu(self):
        try:
            return u.Quantity(self.virt_packet_nus, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_nu:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_energy(self):
        try:
            return u.Quantity(self.virt_packet_energies, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_energy:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return self.virtual_packet_energy / self.time_of_simulation
        except TypeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_luminosity:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def montecarlo_virtual_luminosity(self):
        return (
            self._montecarlo_virtual_luminosity[:-1]
            / self.time_of_simulation.value
        )
