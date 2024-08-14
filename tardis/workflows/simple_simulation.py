import logging
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.io.atom_data.base import AtomData
from tardis.model import SimulationState
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.plasma.assembly.legacy_assembly import assemble_plasma
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral import FormalIntegrator
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.base import MonteCarloTransportSolver
from tardis.util.base import is_notebook
from tardis.workflows.workflow_logging import WorkflowLogging
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.macro_atom.macroatom_solver import MacroAtomSolver

# logging support
logger = logging.getLogger(__name__)


class SimpleSimulation(WorkflowLogging):
    show_progress_bars = is_notebook()
    enable_virtual_packet_logging = False
    log_level = None
    specific_log_level = None

    def __init__(self, configuration):
        super().__init__(configuration, self.log_level, self.specific_log_level)
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
            enable_virtual_packet_logging=self.enable_virtual_packet_logging,
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

        self.opacity_solver = OpacitySolver(
            line_interaction_type=configuration.plasma.line_interaction_type,
            disable_line_scattering=False,
        )
        if configuration.plasma.line_interaction_type == "scatter":
            self.macro_atom_solver = None
        else:
            self.macro_atom_solver = MacroAtomSolver()

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

        luminosity_ratios = (
            (emitted_luminosity / self.luminosity_requested).to(1).value
        )

        estimated_t_inner = (
            self.simulation_state.t_inner
            * luminosity_ratios
            ** self.convergence_strategy.t_inner_update_exponent
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
            return self.consecutive_converges_count >= hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

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

        Returns
        -------
        next_values : dict
            The next values assigned to the simulation state
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

        return next_values

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
            temperature=self.simulation_state.radiation_field_state.temperature,
            dilution_factor=self.simulation_state.radiation_field_state.dilution_factor,
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

        opacity_state = self.opacity_solver.solve(self.plasma_solver)

        if self.macro_atom_solver is None:
            macro_atom_state = None
        else:
            macro_atom_state = self.macro_atom_solver.solve(
                self.plasma_solver,
                self.plasma_solver.atomic_data,
                opacity_state.tau_sobolev,
                self.plasma_solver.stimulated_emission_factor,
            )
        transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
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

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        converged = False
        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                f"\n\tStarting iteration {(self.completed_iterations + 1):d} of {self.total_iterations:d}"
            )
            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates(transport_state)

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(estimated_radfield_properties)

            converged = self.check_convergence(estimated_values)
            self.completed_iterations += 1

            if converged and self.convergence_strategy.stop_if_converged:
                break

        if converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )
        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count, self.virtual_packet_count
        )

        self.initialize_spectrum_solver(
            transport_state,
            virtual_packet_energies,
        )
