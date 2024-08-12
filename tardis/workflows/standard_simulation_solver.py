import logging
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.io.atom_data.base import AtomData
from tardis.io.util import HDFWriterMixin
from tardis.model import SimulationState
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.simulation.base import PlasmaStateStorerMixin
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral import FormalIntegrator
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.base import MonteCarloTransportSolver
from tardis.util.base import is_notebook
from tardis.visualization import ConvergencePlots
from tardis.workflows.workflow_logging import WorkflowLogging

# logging support
logger = logging.getLogger(__name__)


class StandardSimulationSolver(
    WorkflowLogging, PlasmaStateStorerMixin, HDFWriterMixin
):
    hdf_properties = [
        "simulation_state",
        "plasma_solver",
        "transport_solver",
        "iterations_w",
        "iterations_t_rad",
        "iterations_electron_densities",
        "iterations_t_inner",
        "spectrum_solver",
    ]

    def __init__(
        self,
        configuration,
        enable_virtual_packet_logging=False,
        log_level=None,
        specific_log_level=None,
        show_progress_bars=False,
        show_convergence_plots=False,
        convergence_plots_kwargs={},
    ):
        # set up logging
        WorkflowLogging.__init__(
            self,
            configuration=configuration,
            log_level=log_level,
            specific_log_level=specific_log_level,
        )

        self.show_progress_bars = show_progress_bars

        atom_data = self._get_atom_data(configuration)

        # set up states and solvers
        self.simulation_state = SimulationState.from_config(
            configuration,
            atom_data=atom_data,
        )

        self.plasma_solver = assemble_plasma(
            configuration,
            self.simulation_state,
            atom_data=atom_data,
        )

        self.transport_solver = MonteCarloTransportSolver.from_config(
            configuration,
            packet_source=self.simulation_state.packet_source,
            enable_virtual_packet_logging=enable_virtual_packet_logging,
        )

        # Luminosity filter frequencies
        self.luminosity_nu_start = (
            configuration.supernova.luminosity_wavelength_end.to(
                u.Hz, u.spectral()
            )
        )

        if u.isclose(
            configuration.supernova.luminosity_wavelength_start, 0 * u.angstrom
        ):
            self.luminosity_nu_end = np.inf * u.Hz
        else:
            self.luminosity_nu_end = (
                const.c / configuration.supernova.luminosity_wavelength_start
            ).to(u.Hz)

        # montecarlo settings
        self.total_iterations = int(configuration.montecarlo.iterations)

        self.real_packet_count = int(configuration.montecarlo.no_of_packets)

        final_iteration_packet_count = (
            configuration.montecarlo.last_no_of_packets
        )

        if (
            final_iteration_packet_count is None
            or final_iteration_packet_count < 0
        ):
            final_iteration_packet_count = self.real_packet_count

        self.final_iteration_packet_count = int(final_iteration_packet_count)

        self.virtual_packet_count = int(
            configuration.montecarlo.no_of_virtual_packets
        )

        # set up plasma storage
        PlasmaStateStorerMixin.__init__(
            self,
            iterations=self.total_iterations,
            no_of_shells=self.simulation_state.no_of_shells,
        )

        # spectrum settings
        self.integrated_spectrum_settings = configuration.spectrum.integrated
        self.spectrum_solver = SpectrumSolver.from_config(configuration)

        # Convergence settings
        self.consecutive_converges_count = 0
        self.converged = False
        self.completed_iterations = 0
        self.luminosity_requested = (
            configuration.supernova.luminosity_requested.cgs
        )

        # Convergence solvers
        self.convergence_strategy = (
            configuration.montecarlo.convergence_strategy
        )

        self.convergence_solvers = {}
        self.convergence_solvers["t_radiative"] = ConvergenceSolver(
            self.convergence_strategy.t_rad
        )
        self.convergence_solvers["dilution_factor"] = ConvergenceSolver(
            self.convergence_strategy.w
        )
        self.convergence_solvers["t_inner"] = ConvergenceSolver(
            self.convergence_strategy.t_inner
        )

        # Convergence plots
        if show_convergence_plots:
            if not is_notebook():
                raise RuntimeError(
                    "Convergence Plots cannot be displayed in command-line. Set show_convergence_plots "
                    "to False."
                )

            self.convergence_plots = ConvergencePlots(
                iterations=self.total_iterations, **convergence_plots_kwargs
            )
        else:
            self.convergence_plots = None

        if "export_convergence_plots" in convergence_plots_kwargs:
            if not isinstance(
                convergence_plots_kwargs["export_convergence_plots"],
                bool,
            ):
                raise TypeError(
                    "Expected bool in export_convergence_plots argument"
                )
            self.export_convergence_plots = convergence_plots_kwargs[
                "export_convergence_plots"
            ]
        else:
            self.export_convergence_plots = False

    def _get_atom_data(self, configuration):
        """Process atomic data from the configuration

        Parameters
        ----------
        configuration : Configuration
            TARDIS configuration object

        Returns
        -------
        AtomData
            Atomic data object

        Raises
        ------
        ValueError
            If atom data is missing from the configuration
        """
        if "atom_data" in configuration:
            if Path(configuration.atom_data).is_absolute():
                atom_data_fname = Path(configuration.atom_data)
            else:
                atom_data_fname = (
                    Path(configuration.config_dirname) / configuration.atom_data
                )

        else:
            raise ValueError("No atom_data option found in the configuration.")

        logger.info(f"\n\tReading Atomic Data from {atom_data_fname}")

        try:
            atom_data = AtomData.from_hdf(atom_data_fname)
        except TypeError:
            logger.exception(
                "TypeError might be from the use of an old-format of the atomic database, \n"
                "please see https://github.com/tardis-sn/tardis-refdata/tree/master/atom_data"
                " for the most recent version.",
            )
            raise

        return atom_data

    def get_convergence_estimates(self, transport_state):
        """Compute convergence estimates from the transport state

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            Transport state object to compute estimates

        Returns
        -------
        dict
            Convergence estimates
        EstimatedRadiationFieldProperties
            Dilute radiation file and j_blues dataclass
        """
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                transport_state.radfield_mc_estimators,
                transport_state.time_explosion,
                transport_state.time_of_simulation,
                transport_state.geometry_state.volume,
                transport_state.opacity_state.line_list_nu,
            )
        )

        estimated_t_radiative = (
            estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        )
        estimated_dilution_factor = (
            estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor
        )

        emitted_luminosity = calculate_filtered_luminosity(
            transport_state.emitted_packet_nu,
            transport_state.emitted_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )
        absorbed_luminosity = calculate_filtered_luminosity(
            transport_state.reabsorbed_packet_nu,
            transport_state.reabsorbed_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )

        if self.convergence_plots is not None:
            plot_data = {
                "t_inner": [self.simulation_state.t_inner.value, "value"],
                "t_rad": [self.simulation_state.t_radiative, "iterable"],
                "w": [self.simulation_state.dilution_factor, "iterable"],
                "velocity": [self.simulation_state.velocity, "iterable"],
                "Emitted": [emitted_luminosity.value, "value"],
                "Absorbed": [absorbed_luminosity.value, "value"],
                "Requested": [self.luminosity_requested.value, "value"],
            }
            self.update_convergence_plot_data(plot_data)

        self.log_iteration_results(
            emitted_luminosity, absorbed_luminosity, self.luminosity_requested
        )

        luminosity_ratios = (
            (emitted_luminosity / self.luminosity_requested).to(1).value
        )

        estimated_t_inner = (
            self.simulation_state.t_inner
            * luminosity_ratios
            ** self.convergence_strategy.t_inner_update_exponent
        )

        self.log_plasma_state(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_radiative,
            estimated_dilution_factor,
            estimated_t_inner,
        )

        return {
            "t_radiative": estimated_t_radiative,
            "dilution_factor": estimated_dilution_factor,
            "t_inner": estimated_t_inner,
        }, estimated_radfield_properties

    def check_convergence(
        self,
        estimated_values,
    ):
        """Check convergence status for a dict of estimated values

        Parameters
        ----------
        estimated_values : dict
            Estimates to check convergence

        Returns
        -------
        bool
            If convergence has occurred
        """
        convergence_statuses = []

        for key, solver in self.convergence_solvers.items():
            current_value = getattr(self.simulation_state, key)
            estimated_value = estimated_values[key]
            no_of_shells = (
                self.simulation_state.no_of_shells if key != "t_inner" else 1
            )
            convergence_statuses.append(
                solver.get_convergence_status(
                    current_value, estimated_value, no_of_shells
                )
            )

        if np.all(convergence_statuses):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            return self.consecutive_converges_count == hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def update_convergence_plot_data(self, plot_data_dict):
        """Updates convergence plotting data

        Parameters
        ----------
        plot_data_dict : dict
            Dictionary of data to update of the form {"name": [value, item_type]}
        """
        for name, (value, item_type) in plot_data_dict.items():
            self.convergence_plots.fetch_data(
                name=name,
                value=value,
                item_type=item_type,
            )

    def solve_simulation_state(
        self,
        estimated_values,
    ):
        """Update the simulation state with new inputs computed from previous
        iteration estimates.

        Parameters
        ----------
        estimated_values : dict
            Estimated from the previous iterations
        """
        next_values = {}

        for key, solver in self.convergence_solvers.items():
            if (
                key == "t_inner"
                and (self.completed_iterations + 1)
                % self.convergence_strategy.lock_t_inner_cycles
                != 0
            ):
                next_values[key] = getattr(self.simulation_state, key)
            else:
                next_values[key] = solver.converge(
                    getattr(self.simulation_state, key), estimated_values[key]
                )

        self.simulation_state.t_radiative = next_values["t_radiative"]
        self.simulation_state.dilution_factor = next_values["dilution_factor"]
        self.simulation_state.blackbody_packet_source.temperature = next_values[
            "t_inner"
        ]

    def solve_plasma(self, estimated_radfield_properties):
        """Update the plasma solution with the new radiation field estimates

        Parameters
        ----------
        estimated_radfield_properties : EstimatedRadiationFieldProperties
            The radiation field properties to use for updating the plasma

        Raises
        ------
        ValueError
            If the plasma solver radiative rates type is unknown
        """
        radiation_field = DilutePlanckianRadiationField(
            temperature=self.simulation_state.t_radiative,
            dilution_factor=self.simulation_state.dilution_factor,
        )
        update_properties = dict(
            dilute_planckian_radiation_field=radiation_field
        )
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.
        if (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "blackbody"
        ):
            planckian_radiation_field = (
                radiation_field.to_planckian_radiation_field()
            )
            j_blues = planckian_radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma_solver.atomic_data.lines.index
            )
        elif (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "dilute-blackbody"
        ):
            j_blues = radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma_solver.atomic_data.lines.index
            )
        elif (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "detailed"
        ):
            update_properties["j_blues"] = pd.DataFrame(
                estimated_radfield_properties.j_blues,
                index=self.plasma_solver.atomic_data.lines.index,
            )
        else:
            raise ValueError(
                f"radiative_rates_type type unknown - {self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE}"
            )

        self.plasma_solver.update(**update_properties)

    def solve_montecarlo(self, no_of_real_packets, no_of_virtual_packets=0):
        """Solve the MonteCarlo process

        Parameters
        ----------
        no_of_real_packets : int
            Number of real packets to simulate
        no_of_virtual_packets : int, optional
            Number of virtual packets to simulate per interaction, by default 0

        Returns
        -------
        MonteCarloTransportState
            The new transport state after simulation
        ndarray
            Array of unnormalized virtual packet energies in each frequency bin
        """
        transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport_solver.run(
            transport_state,
            iteration=self.completed_iterations,
            total_iterations=self.total_iterations,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        return transport_state, virtual_packet_energies

    def initialize_spectrum_solver(
        self,
        transport_state,
        virtual_packet_energies=None,
    ):
        """Set up the spectrum solver

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            The transport state to init with
        virtual_packet_energies : ndarray, optional
            Array of virtual packet energies binned by frequency, by default None
        """
        # Set up spectrum solver
        self.spectrum_solver.transport_state = transport_state

        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[
                :
            ] = virtual_packet_energies

        if self.integrated_spectrum_settings is not None:
            # Set up spectrum solver integrator
            self.spectrum_solver.integrator_settings = (
                self.integrated_spectrum_settings
            )
            self.spectrum_solver._integrator = FormalIntegrator(
                self.simulation_state, self.plasma_solver, self.transport_solver
            )

    def solve(self):
        """Solve the TARDIS simulation until convergence is reached"""
        converged = False
        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                f"\n\tStarting iteration {(self.completed_iterations + 1):d} of {self.total_iterations:d}"
            )
            self.store_plasma_state(
                self.completed_iterations,
                self.simulation_state.dilution_factor,
                self.simulation_state.t_radiative,
                self.plasma_solver.electron_densities,
                self.simulation_state.t_inner,
            )
            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates(transport_state)

            if self.convergence_plots is not None:
                self.convergence_plots.update()

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(estimated_radfield_properties)

            converged = self.check_convergence(estimated_values)
            self.completed_iterations += 1

            if converged and self.convergence_strategy.stop_if_converged:
                break

        logger.info(f"\n\tStarting final iteration")
        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count, self.virtual_packet_count
        )
        self.store_plasma_state(
            self.completed_iterations,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_radiative,
            self.plasma_solver.electron_densities,
            self.simulation_state.t_inner,
        )
        self.reshape_plasma_state_store(self.completed_iterations)
        if self.convergence_plots is not None:
            self.get_convergence_estimates(transport_state)
            self.convergence_plots.update(
                export_convergence_plots=self.export_convergence_plots,
                last=True,
            )
        self.initialize_spectrum_solver(
            transport_state,
            virtual_packet_energies,
        )
