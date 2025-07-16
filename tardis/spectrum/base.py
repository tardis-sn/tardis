import warnings

import numpy as np
from astropy import units as u

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.spectrum.formal_integral.base import IntegrationError
from tardis.spectrum.spectrum import TARDISSpectrum
from tardis.util.base import (
    quantity_linspace,
)


class SpectrumSolver(HDFWriterMixin):
    hdf_properties = [
        "montecarlo_virtual_luminosity",
        "spectrum_real_packets",
        "spectrum_virtual_packets",
        "spectrum_real_packets_reabsorbed",
        "spectrum_integrated",
    ]

    hdf_name = "spectrum"

    def __init__(
        self, transport_state, spectrum_frequency_grid, integrator_settings
    ):
        self.transport_state = transport_state
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency_grid.value), "erg / s"
        )  # should be init with v_packets_energy_hist
        self._integrator = None
        self.integrator_settings = integrator_settings
        self._spectrum_integrated = None

    def setup_optional_spectra(
        self, transport_state, virtual_packet_luminosity=None, integrator=None
    ):
        """Set up the solver to handle real and virtual spectra

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            The transport state to be used to compute the spectra
        v_packets_energy_hist : np.ndarray
            Virtual packets energy histogram, unnormalized
        integrator : FormalIntegratorSolver, optional
            Integrator to compute the integrated spectrum with, by default None
        """
        self.transport_state = transport_state
        self._montecarlo_virtual_luminosity = (
            virtual_packet_luminosity * u.erg / u.s
        )
        self._integrator = integrator

    def setup_virtual_spectrum(
            self, virtual_packet_luminosity=None
    ):
        """Set up for the virtual spectrum

        Parameters
        ----------
        virtual_packet_luminosity : np.ndarray, optional
            Virtual packet luminosity, unnormalized, by default None
        """
        self._montecarlo_virtual_luminosity = (
            virtual_packet_luminosity * u.erg / u.s
        )
    
    def setup_integrated_spectrum(self, integrator, simulation_state, opacity_state, transport_state, atomic_data, levels):
        """Set up for the integrated spectra

        Parameters
        ----------
        integrator : FormalIntegratorSolver, optional
            Integrator to compute the integrated spectrum with, by default None
        """
        self._integrator = integrator
        self.simulation_state = simulation_state
        self.opacity_state = opacity_state
        self.transport_state = transport_state
        self.atomic_data = atomic_data
        self.levels_index = levels


    @property
    def spectrum_real_packets(self):
        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_emitted_luminosity
        )

    @property
    def spectrum_real_packets_reabsorbed(self):
        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_reabsorbed_luminosity
        )

    @property
    def spectrum_virtual_packets(self):
        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                "SpectrumSolver.spectrum_virtual_packets "
                "is zero. Please run the montecarlo simulation with "
                "no_of_virtual_packets > 0",
                UserWarning,
            )

        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_virtual_luminosity
        )

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None and self.integrator is not None:
            # This was changed from unpacking to specific attributes as compute
            # is not used in calculate_spectrum
            try:
                self._spectrum_integrated = self.integrator.calculate_spectrum(
                    self.spectrum_frequency_grid[:-1],
                    points=self.integrator_settings.points,
                    interpolate_shells=self.integrator_settings.interpolate_shells,
                )
                # self._spectrum_integrated = self.integrator.solve( # TODO: finish
                #     self.simulation_state, 
                #     self.opacity_state, 
                #     self.transport_state, 
                #     self.atomic_data, 
                #     self.levels_index
                # )
            except IntegrationError:
                # if integration is impossible or fails, return an empty spectrum
                warnings.warn(
                    "The FormalIntegrator is not yet implemented for the full "
                    "relativity mode or continuum processes. "
                    "Please run with config option enable_full_relativity: "
                    "False and continuum_processes_enabled: False "
                    "This RETURNS AN EMPTY SPECTRUM!",
                    UserWarning,
                )
                self._spectrum_integrated = TARDISSpectrum(
                    np.array([np.nan, np.nan]) * u.Hz,
                    np.array([np.nan]) * u.erg / u.s,
                )
        return self._spectrum_integrated

    @property
    def integrator(self):
        if self._integrator is None:
            warnings.warn(
                "SpectrumSolver.integrator: "
                "The FormalIntegrator is not yet available."
                "Please run the montecarlo simulation at least once.",
                UserWarning,
            )
        if self.transport_state.enable_full_relativity:
            raise NotImplementedError(
                "The FormalIntegrator is not yet implemented for the full "
                "relativity mode. "
                "Please run with config option enable_full_relativity: "
                "False."
            )
        return self._integrator

    @property
    def montecarlo_reabsorbed_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.transport_state.reabsorbed_packet_nu,
                weights=self.transport_state.reabsorbed_packet_luminosity,
                bins=self.spectrum_frequency_grid,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_emitted_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.transport_state.emitted_packet_nu,
                weights=self.transport_state.emitted_packet_luminosity,
                bins=self.spectrum_frequency_grid,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_virtual_luminosity(self):
        return (
            self._montecarlo_virtual_luminosity[:-1]
            / self.transport_state.time_of_simulation.value
        )

    def solve(self, transport_state):
        """Solve the spectra

        Parameters
        ----------
        transport_state: MonteCarloTransportState
            The transport state to be used to compute the spectra

        Returns
        -------
        tuple(TARDISSpectrum)
            Real, virtual and integrated spectra, if available
        """
        self.transport_state = transport_state

        return (
            self.spectrum_real_packets,
            self.spectrum_virtual_packets,
            self.spectrum_integrated,
        )

    @classmethod
    def from_config(cls, config):
        spectrum_frequency_grid = quantity_linspace(
            config.spectrum.stop.to("Hz", u.spectral()),
            config.spectrum.start.to("Hz", u.spectral()),
            num=config.spectrum.num + 1,
        )

        return cls(
            transport_state=None,
            spectrum_frequency_grid=spectrum_frequency_grid,
            integrator_settings=config.spectrum.integrated,
        )
