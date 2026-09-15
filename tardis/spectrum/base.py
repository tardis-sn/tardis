from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u
from numba.typed import List as TypedList

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.spectrum.formal_integral.base import IntegrationError
from tardis.spectrum.spectrum import TARDISSpectrum
from tardis.spectrum.virtual_packet import calculate_virtual_packet_spectrum
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.modes.nonhomologous.montecarlo_transport_state import (
    MonteCarloTransportStateNonhomologous,
)
from tardis.transport.montecarlo.modes.nonhomologous.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_nonhomologous,
)
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)
from tardis.transport.montecarlo.packets.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_classic,
)
from tardis.util.base import (
    quantity_linspace,
)

if TYPE_CHECKING:
    from tardis.io.configuration.config_reader import ConfigurationNameSpace
    from tardis.model.base import SimulationState
    from tardis.opacities.macro_atom.macroatom_state import (
        LegacyMacroAtomState,
        MacroAtomState,
    )
    from tardis.opacities.opacity_state import OpacityState
    from tardis.plasma.base import BasePlasma
    from tardis.spectrum.formal_integral.formal_integral_solver import (
        FormalIntegralSolver,
    )
    from tardis.transport.montecarlo.modes.classic.solver import (
        MCTransportSolverClassic,
    )
    from tardis.transport.montecarlo.modes.nonhomologous.solver import (
        MCTransportSolverNonhomologous,
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
        self,
        transport_state: (
            MonteCarloTransportState
            | MonteCarloTransportStateNonhomologous
            | None
        ),
        spectrum_frequency_grid: u.Quantity,
        integrator_settings: ConfigurationNameSpace | None,
    ) -> None:
        self.transport_state = transport_state
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency_grid.value), "erg / s"
        )
        self._integrator = None
        self.integrator_settings = integrator_settings
        self._spectrum_integrated = None

    def setup_optional_spectra(
        self,
        transport_state: (
            MonteCarloTransportState | MonteCarloTransportStateNonhomologous
        ),
        integrator: FormalIntegralSolver | None = None,
        simulation_state: SimulationState | None = None,
        transport: (
            MCTransportSolverClassic | MCTransportSolverNonhomologous | None
        ) = None,
        plasma: BasePlasma | None = None,
        opacity_state: OpacityState | None = None,
        macro_atom_state: LegacyMacroAtomState | MacroAtomState | None = None,
    ) -> None:
        """Set up the solver to handle real and virtual spectra

        Parameters
        ----------
        transport_state : MonteCarloTransportState or MonteCarloTransportStateNonhomologous
            Transport state used to compute the spectra.
        integrator : FormalIntegratorSolver, optional
            Integrator used to compute the integrated spectrum.
        simulation_state : SimulationState, optional
            Final simulation state.
        transport : MCTransportSolverClassic or MCTransportSolverNonhomologous, optional
            Completed transport solver used by the formal integral.
        plasma : BasePlasma, optional
            Final plasma state.
        opacity_state : OpacityState, optional
            Final opacity state.
        macro_atom_state : LegacyMacroAtomState or MacroAtomState, optional
            Final macro-atom state.
        """
        self.transport_state = transport_state
        if transport is not None:
            self.solve_virtual_packet_spectrum(
                transport_state,
                transport.montecarlo_configuration,
                transport.rpacket_tracker,
            )
            transport.rpacket_tracker = None
        self._integrator = integrator
        self.simulation_state = simulation_state
        self.opacity_state = opacity_state
        self.transport = transport
        self.plasma = plasma
        self.macro_atom_state = macro_atom_state

    def solve_virtual_packet_spectrum(
        self,
        transport_state: (
            MonteCarloTransportState | MonteCarloTransportStateNonhomologous
        ),
        montecarlo_configuration: MonteCarloConfiguration,
        rpacket_tracker: TypedList[TrackerFull] | None,
    ) -> None:
        """Compute the virtual spectrum from retained real-packet histories.

        Parameters
        ----------
        transport_state : MonteCarloTransportState or MonteCarloTransportStateNonhomologous
            Completed classic or non-homologous transport state.
        montecarlo_configuration : MonteCarloConfiguration
            Configuration used for the completed transport calculation.
        rpacket_tracker : numba.typed.List, optional
            Retained full real-packet histories. Required when virtual packets
            are configured.
        """
        number_of_vpackets = montecarlo_configuration.NUMBER_OF_VPACKETS

        self._montecarlo_virtual_luminosity.value[:] = 0.0
        transport_state.virt_logging = (
            montecarlo_configuration.ENABLE_VPACKET_TRACKING
        )
        if number_of_vpackets <= 0:
            return

        if isinstance(transport_state, MonteCarloTransportState):
            time_explosion = transport_state.time_explosion.cgs.value
            trace_vpacket_volley = trace_vpacket_volley_classic
        elif isinstance(
            transport_state, MonteCarloTransportStateNonhomologous
        ):
            time_explosion = 0.0
            trace_vpacket_volley = trace_vpacket_volley_nonhomologous
        else:
            raise TypeError(
                "Virtual-packet spectra are supported only for classic and "
                "non-homologous Monte Carlo transport states."
            )

        if rpacket_tracker is None:
            raise RuntimeError(
                "Virtual-packet post-processing requires completed "
                "real-packet histories."
            )

        virtual_packet_luminosity, vpacket_tracker = (
            calculate_virtual_packet_spectrum(
                transport_state.packet_collection,
                transport_state.geometry_state_numba,
                time_explosion,
                transport_state.opacity_state_numba,
                montecarlo_configuration,
                self.spectrum_frequency_grid.value,
                rpacket_tracker,
                number_of_vpackets,
                trace_vpacket_volley,
            )
        )
        self._montecarlo_virtual_luminosity.value[:] = virtual_packet_luminosity
        if montecarlo_configuration.ENABLE_VPACKET_TRACKING:
            transport_state.vpacket_tracker = vpacket_tracker

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
                self._spectrum_integrated = self.integrator.solve(
                    self.spectrum_frequency_grid[:-1],
                    self.simulation_state,
                    self.transport,
                    self.opacity_state,
                    self.plasma.atomic_data,
                    self.plasma.electron_densities,
                    self.macro_atom_state,
                )
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
