import logging

import numpy as np
from astropy import units as u
from numba import cuda, set_num_threads

from tardis import constants as const
from tardis.io.logger import montecarlo_tracking as mc_tracker
from tardis.io.util import HDFWriterMixin
from tardis.montecarlo.estimators.radfield_mc_estimators import (
    initialize_estimator_statistics,
)
from tardis.montecarlo.montecarlo_configuration import (
    MonteCarloConfiguration,
    configuration_initialize,
)
from tardis.montecarlo.montecarlo_numba import (
    montecarlo_main_loop,
    numba_config,
)
from tardis.montecarlo.montecarlo_numba.formal_integral import FormalIntegrator
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    NumbaModel,
    opacity_state_initialize,
)
from tardis.montecarlo.montecarlo_numba.r_packet import (
    rpacket_trackers_to_dataframe,
)
from tardis.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
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
        spectrum_frequency,
        virtual_spectrum_spawn_range,
        disable_electron_scattering,
        enable_reflective_inner_boundary,
        enable_full_relativity,
        inner_boundary_albedo,
        line_interaction_type,
        integrator_settings,
        v_packet_settings,
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
        # inject different packets
        self.disable_electron_scattering = disable_electron_scattering
        self.spectrum_frequency = spectrum_frequency
        self.virtual_spectrum_spawn_range = virtual_spectrum_spawn_range
        self.enable_reflective_inner_boundary = enable_reflective_inner_boundary
        self.inner_boundary_albedo = inner_boundary_albedo
        self.enable_full_relativity = enable_full_relativity
        self.line_interaction_type = line_interaction_type
        self.integrator_settings = integrator_settings
        self.v_packet_settings = v_packet_settings
        self.spectrum_method = spectrum_method
        self._integrator = None

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

        # if self.spectrum_method == "integrated":
        #    self.optional_hdf_properties.append("spectrum_integrated")

    def initialize_transport_state(
        self,
        simulation_state,
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
        estimators = initialize_estimator_statistics(
            plasma.tau_sobolevs.shape, gamma_shape
        )

        geometry_state = simulation_state.geometry.to_numba()
        opacity_state = opacity_state_initialize(
            plasma,
            self.line_interaction_type,
            self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            self.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
        )
        transport_state = MonteCarloTransportState(
            packet_collection,
            estimators,
            spectrum_frequency=self.spectrum_frequency,
            geometry_state=geometry_state,
            opacity_state=opacity_state,
        )

        transport_state.enable_full_relativity = self.enable_full_relativity
        transport_state.integrator_settings = self.integrator_settings
        transport_state._integrator = FormalIntegrator(
            simulation_state, plasma, self
        )
        configuration_initialize(
            self.montecarlo_configuration, self, no_of_virtual_packets
        )

        return transport_state

    def run(
        self,
        transport_state,
        time_explosion,
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

        numba_model = NumbaModel(time_explosion.to("s").value)

        number_of_vpackets = self.montecarlo_configuration.NUMBER_OF_VPACKETS

        (
            v_packets_energy_hist,
            last_interaction_tracker,
            vpacket_tracker,
            rpacket_trackers,
        ) = montecarlo_main_loop(
            transport_state.packet_collection,
            transport_state.geometry_state,
            numba_model,
            transport_state.opacity_state,
            self.montecarlo_configuration,
            transport_state.radfield_mc_estimators,
            transport_state.spectrum_frequency.value,
            number_of_vpackets,
            iteration=iteration,
            show_progress_bars=show_progress_bars,
            total_iterations=total_iterations,
            enable_virtual_packet_logging=self.enable_vpacket_tracking,
        )

        transport_state._montecarlo_virtual_luminosity.value[
            :
        ] = v_packets_energy_hist
        transport_state.last_interaction_type = last_interaction_tracker.types
        transport_state.last_interaction_in_nu = last_interaction_tracker.in_nus
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
        # Condition for Checking if RPacket Tracking is enabled
        if self.montecarlo_configuration.ENABLE_RPACKET_TRACKING:
            transport_state.rpacket_tracker = rpacket_trackers

        if self.transport_state.rpacket_tracker is not None:
            self.transport_state.rpacket_tracker_df = (
                rpacket_trackers_to_dataframe(
                    self.transport_state.rpacket_tracker
                )
            )
        transport_state.virt_logging = (
            self.montecarlo_configuration.ENABLE_VPACKET_TRACKING
        )

    def legacy_return(self):
        return (
            self.transport_state.packet_collection.output_nus,
            self.transport_state.packet_collection.output_energies,
            self.transport_state.j_estimator,
            self.transport_state.nu_bar_estimator,
            self.transport_state.last_line_interaction_in_id,
            self.transport_state.last_line_interaction_out_id,
            self.transport_state.last_interaction_type,
            self.transport_state.last_line_interaction_shell_id,
        )

    def get_line_interaction_id(self, line_interaction_type):
        return ["scatter", "downbranch", "macroatom"].index(
            line_interaction_type
        )

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
            numba_config.SIGMA_THOMSON = 1e-200
            # mc_config_module.disable_electron_scattering = True
        else:
            logger.debug("Electron scattering switched on")
            numba_config.SIGMA_THOMSON = const.sigma_T.to("cm^2").value
            # mc_config_module.disable_electron_scattering = False

        spectrum_frequency = quantity_linspace(
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

        montecarlo_configuration.INITIAL_TRACKING_ARRAY_LENGTH = (
            config.montecarlo.tracking.initial_array_length
        )

        return cls(
            spectrum_frequency=spectrum_frequency,
            virtual_spectrum_spawn_range=config.montecarlo.virtual_spectrum_spawn_range,
            enable_reflective_inner_boundary=config.montecarlo.enable_reflective_inner_boundary,
            inner_boundary_albedo=config.montecarlo.inner_boundary_albedo,
            enable_full_relativity=config.montecarlo.enable_full_relativity,
            line_interaction_type=config.plasma.line_interaction_type,
            integrator_settings=config.spectrum.integrated,
            v_packet_settings=config.spectrum.virtual,
            spectrum_method=config.spectrum.method,
            disable_electron_scattering=config.plasma.disable_electron_scattering,
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
