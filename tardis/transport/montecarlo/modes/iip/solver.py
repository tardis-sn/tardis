import logging

from astropy import units as u
from numba import cuda, set_num_threads

import tardis.transport.montecarlo.configuration.constants as constants
from tardis import constants as const
from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.io.logger import montecarlo_tracking as mc_tracker
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
    configuration_initialize,
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


# TODO: refactor this into more parts
class MCTransportSolverIIP(HDFWriterMixin):
    """
    This class modifies the MonteCarloTransportState to solve the radiative
    transfer problem.
    """

    hdf_properties = ["transport_state"]

    hdf_name = "transport"

    def __init__(
        self,
        radfield_prop_solver,
        spectrum_frequency_grid,
        virtual_spectrum_spawn_range,
        enable_full_relativity,
        line_interaction_type,
        spectrum_method,
        packet_source,
        enable_virtual_packet_logging=False,
        enable_rpacket_tracking=False,
        nthreads=1,
        debug_packets=False,
        logger_buffer=1,
        use_gpu=False,
        montecarlo_configuration=None,
    ):
        self.radfield_prop_solver = radfield_prop_solver
        # inject different packets
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self.virtual_spectrum_spawn_range = virtual_spectrum_spawn_range
        self.enable_full_relativity = enable_full_relativity
        self.line_interaction_type = line_interaction_type
        self.spectrum_method = spectrum_method

        self.use_gpu = use_gpu

        self.enable_vpacket_tracking = enable_virtual_packet_logging
        self.enable_rpacket_tracking = enable_rpacket_tracking
        self.montecarlo_configuration = montecarlo_configuration

        self.packet_source = packet_source

        # Setting up the Tracking array for storing all the RPacketTracker instances
        self.rpacket_tracker = None

        # Set number of threads
        self.nthreads = nthreads

        # set up logger based on config
        mc_tracker.DEBUG_MODE = debug_packets
        mc_tracker.BUFFER = logger_buffer

    def initialize_transport_state(
        self,
        simulation_state,
        opacity_state,
        macro_atom_state,
        plasma,
        no_of_packets,
        no_of_virtual_packets=0,
        iteration=0,
    ):
        if not plasma.continuum_interaction_species.empty:
            if plasma.gamma is not None:
                n_levels_bf_species_by_n_cells_tuple = plasma.gamma.shape
            else:
                n_levels_bf_species_by_n_cells_tuple = plasma.phi_lucy.shape
        else:
            n_levels_bf_species_by_n_cells_tuple = (0, 0)

        packet_collection = self.packet_source.create_packets(
            no_of_packets, seed_offset=iteration
        )

        # IIP mode: continuum processes always enabled
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True

        geometry_state = simulation_state.geometry.to_numba()
        opacity_state_numba = opacity_state.to_numba(
            macro_atom_state,
            self.line_interaction_type,
        )
        opacity_state_numba = opacity_state_numba[
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
        ]

        transport_state = MonteCarloTransportState(
            packet_collection,
            geometry_state=geometry_state,
            opacity_state=opacity_state_numba,
            time_explosion=simulation_state.time_explosion,
            n_levels_bf_species_by_n_cells_tuple=n_levels_bf_species_by_n_cells_tuple,
        )

        # IIP mode: full relativity always enabled
        transport_state.enable_full_relativity = True

        configuration_initialize(
            self.montecarlo_configuration, self, no_of_virtual_packets
        )

        return transport_state

    def run(
        self,
        transport_state,
        show_progress_bars=True,
    ):
        """
        Run the Monte Carlo calculation using IIP mode (continuum always enabled).

        Parameters
        ----------
        transport_state : tardis.transport.montecarlo.transport_state.TransportState
            Transport state containing all the data needed for the Monte Carlo simulation
        show_progress_bars : bool
            Show progress bars

        Returns
        -------
        v_packets_energy_hist : ndarray
            Histogram of energy from virtual packets
        """
        set_num_threads(self.nthreads)
        self.transport_state = transport_state

        number_of_vpackets = self.montecarlo_configuration.NUMBER_OF_VPACKETS
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

        # IIP mode: returns 3 estimator objects (bulk, line, continuum)
        (
            estimators_bulk,
            estimators_line,
            estimators_continuum,
        ) = montecarlo_transport(
            transport_state.packet_collection,
            transport_state.geometry_state,
            transport_state.time_explosion.cgs.value,
            transport_state.opacity_state,
            self.montecarlo_configuration,
            transport_state.n_levels_bf_species_by_n_cells_tuple,
            trackers_list,
            show_progress_bars=show_progress_bars,
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

        # IIP mode does not currently track virtual packets in montecarlo_transport

    @classmethod
    def from_config(
        cls, config, packet_source, enable_virtual_packet_logging=False
    ):
        """
        Create a new MontecarloTransport instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        virtual_packet_logging : bool

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
            virtual_spectrum_spawn_range=config.montecarlo.virtual_spectrum_spawn_range,
            enable_full_relativity=config.montecarlo.enable_full_relativity,
            line_interaction_type=config.plasma.line_interaction_type,
            spectrum_method=config.spectrum.method,
            packet_source=packet_source,
            debug_packets=config.montecarlo.debug_packets,
            logger_buffer=config.montecarlo.logger_buffer,
            enable_virtual_packet_logging=(
                config.spectrum.virtual.virtual_packet_logging
                | enable_virtual_packet_logging
            ),
            enable_rpacket_tracking=config.montecarlo.tracking.track_rpacket,
            nthreads=config.montecarlo.nthreads,
            use_gpu=use_gpu,
            montecarlo_configuration=montecarlo_configuration,
        )
