import logging
from dataclasses import dataclass
from itertools import pairwise
from typing import Self

import numpy as np
import numpy.typing as npt
from astropy import units as u
from numba import cuda, set_num_threads
from numba.typed import List

import tardis.transport.montecarlo.configuration.constants as constants
from tardis import constants as const
from tardis.io.configuration.config_reader import Configuration
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.io.logger import montecarlo_tracking as mc_tracker
from tardis.opacities.continuum.continuum_state import ContinuumOpacityState
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
    configuration_initialize,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
    init_estimators_continuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
    init_estimators_line,
)
from tardis.transport.montecarlo.estimators.mc_rad_field_solver import (
    MCRadiationFieldPropertiesSolver,
)
from tardis.transport.montecarlo.modes.iip.montecarlo_transport import (
    montecarlo_transport,
)
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.packet_source.base import BasePacketSource
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full_util import (
    generate_tracker_full_list,
    tracker_full_df2tracker_last_interaction_df,
    trackers_full_to_df,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction_util import (
    generate_tracker_last_interaction_list,
    trackers_last_interaction_to_df,
)
from tardis.transport.montecarlo.progress_bars import (
    refresh_packet_pbar,
    reset_packet_pbar,
    update_iterations_pbar,
)
from tardis.util.base import (
    quantity_linspace,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class EstimatorBatches:
    """Estimators of independent, consecutive packet ranges.

    Every packet has its own random-number seed, so each batch is an
    independent Monte Carlo estimate, and the batches sum to the totals.
    """

    estimators_bulk: tuple[EstimatorsBulk, ...]
    estimators_continuum: tuple[EstimatorsContinuum, ...]
    mean_intensity_blueward: tuple[npt.NDArray[np.float64], ...]
    line_idx: npt.NDArray[np.int64]
    packet_counts: npt.NDArray[np.int64]


# TODO: refactor this into more parts
class MCTransportSolverIIP(HDFWriterMixin):
    """Solve Type IIP Monte Carlo radiative transfer.

    Build and evolve the continuum-enabled transport state.
    """

    hdf_properties = ["transport_state"]

    hdf_name = "transport"

    def __init__(
        self,
        radfield_prop_solver: MCRadiationFieldPropertiesSolver,
        spectrum_frequency_grid: u.Quantity,
        enable_full_relativity: bool,
        line_interaction_type: str,
        spectrum_method: str,
        packet_source: BasePacketSource,
        enable_rpacket_tracking: bool = False,
        nthreads: int = 1,
        debug_packets: bool = False,
        logger_buffer: int = 1,
        use_gpu: bool = False,
        montecarlo_configuration: MonteCarloConfiguration | None = None,
    ) -> None:
        self.radfield_prop_solver = radfield_prop_solver
        # inject different packets
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self.enable_full_relativity = enable_full_relativity
        self.line_interaction_type = line_interaction_type
        self.spectrum_method = spectrum_method

        self.use_gpu = use_gpu

        self.enable_rpacket_tracking = enable_rpacket_tracking
        self.montecarlo_configuration = montecarlo_configuration

        self.packet_source = packet_source

        # Set number of threads
        self.nthreads = nthreads

        # set up logger based on config
        mc_tracker.DEBUG_MODE = debug_packets
        mc_tracker.BUFFER = logger_buffer

    def initialize_transport_state(
        self,
        simulation_state: object,
        opacity_state: object,
        macro_atom_state: object | None,
        continuum_state: ContinuumOpacityState,
        no_of_packets: int,
        iteration: int = 0,
    ) -> MonteCarloTransportState:
        """Initialize transport state from explicit opacity inputs."""
        n_levels_bf_species_by_n_cells_tuple = (
            len(continuum_state.level2continuum_idx),
            len(opacity_state.electron_density),
        )

        packet_collection = self.packet_source.create_packets(
            no_of_packets, seed_offset=iteration
        )

        # IIP mode: continuum processes always enabled
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True

        geometry_state = simulation_state.geometry.to_numba()
        opacity_state_numba = opacity_state.to_numba(
            macro_atom_state,
            self.line_interaction_type,
            continuum_processes_enabled=True,
        )
        continuum_state_numba = continuum_state.to_numba(
            opacity_state.t_electrons, macro_atom_state
        )

        transport_state = MonteCarloTransportState(
            packet_collection,
            geometry_state_numba=geometry_state,
            opacity_state_numba=opacity_state_numba,
            time_explosion=simulation_state.time_explosion,
            n_levels_bf_species_by_n_cells_tuple=n_levels_bf_species_by_n_cells_tuple,
        )
        transport_state.continuum_state_numba = continuum_state_numba

        # IIP mode: full relativity always enabled
        transport_state.enable_full_relativity = True

        configuration_initialize(self.montecarlo_configuration, self)

        return transport_state

    def transport_packet_batches(
        self,
        transport_state: MonteCarloTransportState,
        trackers_list: List,
        show_progress_bars: bool,
        estimator_batch_count: int,
        batch_line_idx: npt.NDArray[np.int64] | None,
    ) -> tuple[EstimatorsBulk, EstimatorsLine, EstimatorsContinuum]:
        """Transport consecutive packet ranges and sum their estimators.

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            Transport state whose packet collection receives the outputs.
        trackers_list : numba.typed.List
            Per-packet trackers in packet order.
        show_progress_bars : bool
            Show progress bars.
        estimator_batch_count : int
            Number of packet ranges. Only one range keeps no batches.
        batch_line_idx : numpy.ndarray or None
            Line positions whose blueward mean intensity is kept per batch.

        Returns
        -------
        tuple
            Summed bulk, line, and continuum estimators.
        """
        packet_collection = transport_state.packet_collection
        number_of_packets = packet_collection.number_of_packets
        transport_arguments = (
            transport_state.geometry_state_numba,
            transport_state.time_explosion.cgs.value,
            transport_state.opacity_state_numba,
            transport_state.continuum_state_numba,
            self.montecarlo_configuration,
            transport_state.n_levels_bf_species_by_n_cells_tuple,
        )
        transport_state.estimator_batches = None
        if estimator_batch_count <= 1:
            return montecarlo_transport(
                packet_collection,
                *transport_arguments,
                trackers_list,
                show_progress_bars=show_progress_bars,
            )

        n_cells = len(transport_state.geometry_state_numba.r_inner)
        estimators_bulk = init_estimators_bulk(n_cells)
        estimators_line = init_estimators_line(
            transport_state.opacity_state_numba.tau_sobolev.shape
        )
        estimators_continuum = init_estimators_continuum(
            transport_state.n_levels_bf_species_by_n_cells_tuple, n_cells
        )
        batch_edge_idx = np.linspace(
            0, number_of_packets, estimator_batch_count + 1
        ).astype(np.int64)
        batch_bulk, batch_continuum, batch_intensity, batch_counts = (
            [],
            [],
            [],
            [],
        )
        for batch_start_idx, batch_stop_idx in pairwise(batch_edge_idx):
            if batch_stop_idx == batch_start_idx:
                continue
            batch_packets = PacketCollection(
                packet_collection.initial_radii[batch_start_idx:batch_stop_idx],
                packet_collection.initial_nus[batch_start_idx:batch_stop_idx],
                packet_collection.initial_mus[batch_start_idx:batch_stop_idx],
                packet_collection.initial_energies[
                    batch_start_idx:batch_stop_idx
                ],
                packet_collection.packet_seeds[batch_start_idx:batch_stop_idx],
                packet_collection.radiation_field_luminosity,
            )
            bulk, line, continuum = montecarlo_transport(
                batch_packets,
                *transport_arguments,
                trackers_list[batch_start_idx:batch_stop_idx],
                show_progress_bars=show_progress_bars,
            )
            packet_collection.output_nus[batch_start_idx:batch_stop_idx] = (
                batch_packets.output_nus
            )
            packet_collection.output_energies[
                batch_start_idx:batch_stop_idx
            ] = batch_packets.output_energies
            estimators_bulk.increment(bulk)
            estimators_line.increment(line)
            estimators_continuum.increment(continuum)
            batch_bulk.append(bulk)
            batch_continuum.append(continuum)
            batch_intensity.append(
                line.mean_intensity_blueward.copy()
                if batch_line_idx is None
                else line.mean_intensity_blueward[batch_line_idx]
            )
            batch_counts.append(batch_stop_idx - batch_start_idx)

        transport_state.estimator_batches = EstimatorBatches(
            tuple(batch_bulk),
            tuple(batch_continuum),
            tuple(batch_intensity),
            (
                np.arange(estimators_line.mean_intensity_blueward.shape[0])
                if batch_line_idx is None
                else np.asarray(batch_line_idx, dtype=np.int64)
            ),
            np.asarray(batch_counts, dtype=np.int64),
        )
        return estimators_bulk, estimators_line, estimators_continuum

    def run(
        self,
        transport_state: MonteCarloTransportState,
        show_progress_bars: bool = True,
        estimator_batch_count: int = 1,
        batch_line_idx: npt.NDArray[np.int64] | None = None,
    ) -> None:
        """
        Run the Monte Carlo calculation using IIP mode (continuum always enabled).

        Parameters
        ----------
        transport_state : tardis.transport.montecarlo.transport_state.TransportState
            Transport state containing all the data needed for the Monte Carlo simulation
        show_progress_bars : bool
            Show progress bars
        estimator_batch_count : int, optional
            Number of consecutive packet ranges transported separately. With
            more than one, ``transport_state.estimator_batches`` holds each
            range's estimators so their Monte Carlo noise can be measured;
            the summed estimators equal those of a single transport up to the
            order of floating-point summation.
        batch_line_idx : numpy.ndarray, optional
            Line positions whose blueward mean intensity is kept for each
            batch. All lines are kept when omitted.
        """
        set_num_threads(self.nthreads)
        self.transport_state = transport_state

        number_of_rpackets = len(transport_state.packet_collection.initial_nus)

        if self.enable_rpacket_tracking:
            trackers_list = generate_tracker_full_list(
                number_of_rpackets,
                self.montecarlo_configuration.INITIAL_TRACKING_ARRAY_LENGTH,
            )
        else:
            # Initialize the last interaction tracker list directly
            trackers_list = generate_tracker_last_interaction_list(
                number_of_rpackets
            )

        # Reset packet progress bar for this iteration
        if show_progress_bars:
            reset_packet_pbar(number_of_rpackets)

        estimators_bulk, estimators_line, estimators_continuum = (
            self.transport_packet_batches(
                transport_state,
                trackers_list,
                show_progress_bars,
                estimator_batch_count,
                batch_line_idx,
            )
        )

        # Attach estimators to transport state
        transport_state.estimators_bulk = estimators_bulk
        transport_state.estimators_line = estimators_line
        transport_state.estimators_continuum = (
            estimators_continuum  # IIP mode specific
        )

        # Last interaction trackers are already populated directly in the list
        # No finalization needed with direct list approach

        update_iterations_pbar(1)
        refresh_packet_pbar()

        # Need to change the implementation of rpacket_trackers_to_dataframe
        # Such that it also takes of the case of
        # RPacketLastInteractionTracker
        if self.enable_rpacket_tracking:
            self.transport_state.tracker_full_df = trackers_full_to_df(
                trackers_list
            )
            self.transport_state.tracker_last_interaction_df = (
                tracker_full_df2tracker_last_interaction_df(
                    self.transport_state.tracker_full_df
                )
            )
        else:
            self.transport_state.tracker_full_df = None
            self.transport_state.tracker_last_interaction_df = (
                trackers_last_interaction_to_df(trackers_list)
            )

    @classmethod
    def from_config(
        cls, config: Configuration, packet_source: BasePacketSource
    ) -> Self:
        """
        Create a new MontecarloTransport instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        Returns
        -------
        MontecarloTransport
        """
        if config.plasma.disable_electron_scattering:
            logger.warning(
                "Disabling electron scattering - this is not physical."
                "Likely bug in formal integral - "
                "will not give same results."
            )
            constants.SIGMA_THOMSON = 1e-200
        else:
            logger.debug("Electron scattering switched on")
            constants.SIGMA_THOMSON = const.sigma_T.to("cm^2").value

        spectrum_frequency_grid = quantity_linspace(
            config.spectrum.stop.to("Hz", u.spectral()),
            config.spectrum.start.to("Hz", u.spectral()),
            num=config.spectrum.num + 1,
        )
        running_mode = config.spectrum.integrated.compute.upper()

        if running_mode == "GPU":
            if cuda.is_available():
                use_gpu = True
            else:
                raise ValueError(
                    """The GPU option was selected for the formal_integral,
                    but no CUDA GPU is available."""
                )
        elif running_mode == "AUTOMATIC":
            use_gpu = bool(cuda.is_available())
        elif running_mode == "CPU":
            use_gpu = False
        else:
            raise ValueError(
                """An invalid option for compute was passed. The three
                valid values are 'GPU', 'CPU', and 'Automatic'."""
            )

        montecarlo_configuration = MonteCarloConfiguration()

        montecarlo_configuration.DISABLE_LINE_SCATTERING = (
            config.plasma.disable_line_scattering
        )

        montecarlo_configuration.DISABLE_ELECTRON_SCATTERING = (
            config.plasma.disable_electron_scattering
        )

        montecarlo_configuration.INITIAL_TRACKING_ARRAY_LENGTH = (
            config.montecarlo.tracking.initial_array_length
        )

        radfield_prop_solver = MCRadiationFieldPropertiesSolver(
            config.plasma.w_epsilon
        )

        return cls(
            radfield_prop_solver=radfield_prop_solver,
            spectrum_frequency_grid=spectrum_frequency_grid,
            enable_full_relativity=config.montecarlo.enable_full_relativity,
            line_interaction_type=config.plasma.line_interaction_type,
            spectrum_method=config.spectrum.method,
            packet_source=packet_source,
            debug_packets=config.montecarlo.debug_packets,
            logger_buffer=config.montecarlo.logger_buffer,
            enable_rpacket_tracking=config.montecarlo.tracking.track_rpacket,
            nthreads=config.montecarlo.nthreads,
            use_gpu=use_gpu,
            montecarlo_configuration=montecarlo_configuration,
        )
