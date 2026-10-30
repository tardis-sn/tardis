import warnings
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import ConfigurationNameSpace
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.model.base import SimulationState
from tardis.opacities.macro_atom.macroatom_state import (
    LegacyMacroAtomState,
    MacroAtomState,
)
from tardis.opacities.opacity_state import OpacityState
from tardis.plasma.base import BasePlasma
from tardis.spectrum.formal_integral.base import IntegrationError
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.spectrum import TARDISSpectrum
from tardis.spectrum.virtual_packet import (
    calculate_virtual_packet_spectrum,
    create_vpacket_collections,
    tracker_full_df_to_event_arrays,
)
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.modes.classic.solver import (
    MCTransportSolverClassic,
)
from tardis.transport.montecarlo.modes.nonhomologous.montecarlo_transport_state import (
    MonteCarloTransportStateNonhomologous,
)
from tardis.transport.montecarlo.modes.nonhomologous.solver import (
    MCTransportSolverNonhomologous,
)
from tardis.transport.montecarlo.modes.nonhomologous.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_nonhomologous,
)
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
    consolidate_vpacket_tracker,
)
from tardis.transport.montecarlo.packets.trackers import (
    interaction_type_to_categorical,
)
from tardis.transport.montecarlo.packets.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_classic,
)
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
        self,
        transport_state: (
            MonteCarloTransportState
            | MonteCarloTransportStateNonhomologous
            | None
        ),
        spectrum_frequency_grid: u.Quantity,
        integrator_settings: ConfigurationNameSpace | None,
        virtual_packet_logging: bool = False,
        vpacket_spawn_start_frequency: u.Quantity = 0.0 * u.Hz,
        vpacket_spawn_end_frequency: u.Quantity = 1e200 * u.Hz,
        vpacket_tau_russian: float = 10.0,
        vpacket_survival_probability: float = 0.0,
    ) -> None:
        """
        Initialize the spectrum solver.

        Parameters
        ----------
        transport_state : MonteCarloTransportState or MonteCarloTransportStateNonhomologous, optional
            Transport state used to compute the real-packet spectra.
        spectrum_frequency_grid : astropy.units.Quantity
            Frequency-bin edges for the real- and virtual-packet histograms.
        integrator_settings : ConfigurationNameSpace, optional
            Formal integral settings.
        virtual_packet_logging : bool, optional
            Keep the escaping virtual packets in ``vpacket_tracker``.
        vpacket_spawn_start_frequency : astropy.units.Quantity, optional
            Lowest real-packet frequency that spawns virtual packets.
        vpacket_spawn_end_frequency : astropy.units.Quantity, optional
            Highest real-packet frequency that spawns virtual packets.
        vpacket_tau_russian : float, optional
            Optical depth above which virtual packets are Russian-rouletted.
        vpacket_survival_probability : float, optional
            Probability that a Russian-rouletted virtual packet survives.
        """
        self.transport_state = transport_state
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency_grid.value), "erg / s"
        )
        self._integrator = None
        self.integrator_settings = integrator_settings
        self._spectrum_integrated = None
        self.virtual_packet_logging = virtual_packet_logging
        self.vpacket_spawn_start_frequency = vpacket_spawn_start_frequency
        self.vpacket_spawn_end_frequency = vpacket_spawn_end_frequency
        self.vpacket_tau_russian = vpacket_tau_russian
        self.vpacket_survival_probability = vpacket_survival_probability
        self.vpacket_tracker: VPacketCollection | None = None

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
        """Set up the solver to compute the integrated spectrum

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
        # A virtual spectrum from an earlier run is stale for this transport
        self._montecarlo_virtual_luminosity.value[:] = 0.0
        self.vpacket_tracker = None
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
        tracker_full_df: pd.DataFrame | None,
        number_of_vpackets: int,
    ) -> None:
        """Compute the virtual spectrum from full real-packet histories.

        Parameters
        ----------
        transport_state : MonteCarloTransportState or MonteCarloTransportStateNonhomologous
            Completed classic or non-homologous transport state.
        montecarlo_configuration : MonteCarloConfiguration
            Configuration used for the completed transport calculation.
        tracker_full_df : pandas.DataFrame or None
            Full real-packet histories of the completed transport, e.g.
            ``transport_state.tracker_full_df`` or one loaded from disk. Must
            include the ``next_line_id`` column.
        number_of_vpackets : int
            Number of virtual packets spawned per real-packet interaction.
            Must be positive.

        Raises
        ------
        ValueError
            If ``number_of_vpackets`` is not positive, or the full real-packet
            histories are missing or lack a needed column.
        TypeError
            If ``transport_state`` is not a classic or non-homologous state.
        """
        if number_of_vpackets <= 0:
            raise ValueError(
                "Virtual-packet spectra require number_of_vpackets > 0, got "
                f"{number_of_vpackets}."
            )
        if tracker_full_df is None:
            raise ValueError(
                "Virtual-packet post-processing requires full real-packet "
                "histories; run the final iteration with "
                "record_full_histories=True."
            )
        tracker_events = tracker_full_df_to_event_arrays(tracker_full_df)

        if isinstance(transport_state, MonteCarloTransportState):
            time_explosion = transport_state.time_explosion.cgs.value
            trace_vpacket_volley = trace_vpacket_volley_classic
        elif isinstance(transport_state, MonteCarloTransportStateNonhomologous):
            time_explosion = 0.0
            trace_vpacket_volley = trace_vpacket_volley_nonhomologous
        else:
            raise TypeError(
                "Virtual-packet spectra are supported only for classic and "
                "non-homologous Monte Carlo transport states."
            )

        vpacket_spawn_start_frequency = (
            self.vpacket_spawn_start_frequency.to_value(u.Hz)
        )
        vpacket_spawn_end_frequency = self.vpacket_spawn_end_frequency.to_value(
            u.Hz
        )
        vpacket_collections = create_vpacket_collections(
            len(tracker_events.packet_id),
            self.spectrum_frequency_grid.value,
            vpacket_spawn_start_frequency,
            vpacket_spawn_end_frequency,
            number_of_vpackets,
        )
        virtual_packet_luminosity = calculate_virtual_packet_spectrum(
            transport_state.geometry_state_numba,
            time_explosion,
            transport_state.opacity_state_numba,
            montecarlo_configuration,
            self.spectrum_frequency_grid.value,
            tracker_events,
            vpacket_collections,
            self.vpacket_tau_russian,
            self.vpacket_survival_probability,
            trace_vpacket_volley,
        )
        self._montecarlo_virtual_luminosity.value[:] = virtual_packet_luminosity
        if self.virtual_packet_logging:
            self.vpacket_tracker = consolidate_vpacket_tracker(
                vpacket_collections,
                self.spectrum_frequency_grid.value,
                vpacket_spawn_start_frequency,
                vpacket_spawn_end_frequency,
            )

    def vpacket_tracker_to_df(self) -> pd.DataFrame:
        """
        Convert the logged virtual packets to a DataFrame.

        The interaction columns share their names and dtypes with the
        real-packet ``tracker_last_interaction_df``; each row describes the
        real-packet event that spawned the virtual packet.

        Returns
        -------
        pandas.DataFrame
            One row per escaping virtual packet, indexed by ``vpacket_id``,
            with columns ``last_interaction_type``, ``radius``, ``shell_id``,
            ``before_nu``, ``line_absorb_id`` and ``line_emit_id`` for the
            spawning event, and ``nus``, ``energies``, ``initial_mus`` and
            ``initial_rs`` for the virtual packet.

        Raises
        ------
        ValueError
            If no virtual packets were logged.
        """
        if self.vpacket_tracker is None:
            raise ValueError(
                "No virtual packets were logged; enable virtual_packet_logging "
                "and run with no_of_virtual_packets > 0."
            )
        return self._vpacket_tracker_df(
            interaction_type_to_categorical(
                self.vpacket_tracker.last_interaction_type
            )
        )

    def _vpacket_tracker_df(
        self, last_interaction_type: npt.ArrayLike
    ) -> pd.DataFrame:
        """
        Build the logged virtual packet DataFrame.

        Parameters
        ----------
        last_interaction_type : array-like
            Interaction type column, as names or ``InteractionType`` codes.

        Returns
        -------
        pandas.DataFrame
            One row per escaping virtual packet, indexed by ``vpacket_id``.
        """
        vpacket_tracker = self.vpacket_tracker
        return pd.DataFrame(
            {
                "last_interaction_type": last_interaction_type,
                "radius": vpacket_tracker.last_interaction_in_r,
                "shell_id": vpacket_tracker.last_interaction_shell_id,
                "before_nu": vpacket_tracker.last_interaction_in_nu,
                "line_absorb_id": pd.array(
                    vpacket_tracker.last_interaction_in_id, dtype="int64"
                ),
                "line_emit_id": pd.array(
                    vpacket_tracker.last_interaction_out_id, dtype="int64"
                ),
                "nus": vpacket_tracker.nus,
                "energies": vpacket_tracker.energies,
                "initial_mus": vpacket_tracker.initial_mus,
                "initial_rs": vpacket_tracker.initial_rs,
            },
            index=pd.RangeIndex(len(vpacket_tracker.nus), name="vpacket_id"),
        )

    def get_properties(self) -> dict[str, Any]:
        """
        Get properties for HDF storage, adding logged virtual packets.

        Returns
        -------
        dict
            Property names and values; includes ``vpacket_tracker_df`` when
            virtual packets were logged.
        """
        properties = super().get_properties()
        if self.vpacket_tracker is not None:
            # The default fixed HDF format cannot store categoricals, so the
            # interaction type is written as InteractionType codes
            properties["vpacket_tracker_df"] = self._vpacket_tracker_df(
                self.vpacket_tracker.last_interaction_type
            )
        return properties

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

        # The spawn range is given in wavelength, so its end is the lowest
        # spawning frequency
        vpacket_spawn_range = config.montecarlo.virtual_spectrum_spawn_range
        return cls(
            transport_state=None,
            spectrum_frequency_grid=spectrum_frequency_grid,
            integrator_settings=config.spectrum.integrated,
            virtual_packet_logging=config.spectrum.virtual.virtual_packet_logging,
            vpacket_spawn_start_frequency=vpacket_spawn_range.end.to(
                u.Hz, equivalencies=u.spectral()
            ),
            vpacket_spawn_end_frequency=vpacket_spawn_range.start.to(
                u.Hz, equivalencies=u.spectral()
            ),
            vpacket_tau_russian=config.spectrum.virtual.tau_russian,
            vpacket_survival_probability=config.spectrum.virtual.survival_probability,
        )
