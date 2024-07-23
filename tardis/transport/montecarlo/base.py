import logging

from astropy import units as u
from numba import cuda, set_num_threads

import tardis.transport.montecarlo.configuration.constants as constants
from tardis import constants as const
from tardis.io.logger import montecarlo_tracking as mc_tracker
from tardis.io.util import HDFWriterMixin
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
    configuration_initialize,
)
from tardis.transport.montecarlo.estimators.mc_rad_field_solver import (
    MCRadiationFieldPropertiesSolver,
)
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    initialize_estimator_statistics,
)
from tardis.transport.montecarlo.montecarlo_main_loop import (
    montecarlo_main_loop,
)
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.opacities.opacity_state import (
    opacity_state_to_numba,
)
from tardis.transport.montecarlo.packet_trackers import (
    generate_rpacket_tracker_list,
    generate_rpacket_last_interaction_tracker_list,
    rpacket_trackers_to_dataframe,
)
from tardis.util.base import (
    quantity_linspace,
    refresh_packet_pbar,
    update_iterations_pbar,
)

logger = logging.getLogger(__name__)


# TODO: refactor this into more parts
class MonteCarloTransportSolver(HDFWriterMixin):
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
            gamma_shape = plasma.gamma.shape
        else:
            gamma_shape = (0, 0)

        packet_collection = self.packet_source.create_packets(
            no_of_packets, seed_offset=iteration
        )

        geometry_state = simulation_state.geometry.to_numba()

        opacity_state_numba = opacity_state_to_numba(
            opacity_state, macro_atom_state, self.line_interaction_type
        )
        opacity_state_numba = opacity_state_numba[
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
        ]

        estimators = initialize_estimator_statistics(
            opacity_state_numba.tau_sobolev.shape, gamma_shape
        )

        transport_state = MonteCarloTransportState(
            packet_collection,
            estimators,
            geometry_state=geometry_state,
            opacity_state=opacity_state_numba,
            time_explosion=simulation_state.time_explosion,
        )

        transport_state.enable_full_relativity = (
            self.montecarlo_configuration.ENABLE_FULL_RELATIVITY
        )

        configuration_initialize(
            self.montecarlo_configuration, self, no_of_virtual_packets
        )

        return transport_state

    def run(
        self,
        transport_state,
        iteration=0,
        total_iterations=0,
        show_progress_bars=True,
    ):
        """
        Run the montecarlo calculation

        Parameters
        ----------
        model : tardis.model.SimulationState
        plasma : tardis.plasma.BasePlasma
        no_of_packets : int
        no_of_virtual_packets : int
        total_iterations : int
            The total number of iterations in the simulation.

        Returns
        -------
        None
        """
        set_num_threads(self.nthreads)
        self.transport_state = transport_state

        number_of_vpackets = self.montecarlo_configuration.NUMBER_OF_VPACKETS
        number_of_rpackets = len(transport_state.packet_collection.initial_nus)

        if self.enable_rpacket_tracking:
            transport_state.rpacket_tracker = generate_rpacket_tracker_list(
                number_of_rpackets,
                self.montecarlo_configuration.INITIAL_TRACKING_ARRAY_LENGTH,
            )
        else:
            transport_state.rpacket_tracker = (
                generate_rpacket_last_interaction_tracker_list(
                    number_of_rpackets
                )
            )

        (
            v_packets_energy_hist,
            last_interaction_tracker,
            vpacket_tracker,
        ) = montecarlo_main_loop(
            transport_state.packet_collection,
            transport_state.geometry_state,
            transport_state.time_explosion.cgs.value,
            transport_state.opacity_state,
            self.montecarlo_configuration,
            transport_state.radfield_mc_estimators,
            self.spectrum_frequency_grid.value,
            transport_state.rpacket_tracker,
            number_of_vpackets,
            iteration=iteration,
            show_progress_bars=show_progress_bars,
            total_iterations=total_iterations,
        )

        transport_state.last_interaction_type = last_interaction_tracker.types
        transport_state.last_interaction_in_nu = last_interaction_tracker.in_nus
        transport_state.last_interaction_in_r = last_interaction_tracker.in_rs
        transport_state.last_line_interaction_in_id = (
            last_interaction_tracker.in_ids
        )
        transport_state.last_line_interaction_out_id = (
            last_interaction_tracker.out_ids
        )
        transport_state.last_line_interaction_shell_id = (
            last_interaction_tracker.shell_ids
        )

        if self.montecarlo_configuration.ENABLE_VPACKET_TRACKING and (
            number_of_vpackets > 0
        ):
            transport_state.vpacket_tracker = vpacket_tracker

        update_iterations_pbar(1)
        refresh_packet_pbar()

        # Need to change the implementation of rpacket_trackers_to_dataframe
        # Such that it also takes of the case of
        # RPacketLastInteractionTracker
        if self.enable_rpacket_tracking:
            self.transport_state.rpacket_tracker_df = (
                rpacket_trackers_to_dataframe(
                    self.transport_state.rpacket_tracker
                )
            )
        transport_state.virt_logging = (
            self.montecarlo_configuration.ENABLE_VPACKET_TRACKING
        )

        return v_packets_energy_hist

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
