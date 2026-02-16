import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
    ContinuumMacroAtomSolver,
)
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.plasma.assembly import PlasmaSolverFactory
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    MCContinuumPropertiesSolver,
)
from tardis.transport.montecarlo.modes.classic.solver import (
    MCTransportSolverClassic,
)
from tardis.transport.montecarlo.progress_bars import initialize_iterations_pbar
from tardis.util.environment import Environment
from tardis.workflows.workflow_logging import WorkflowLogging

# logging support
logger = logging.getLogger(__name__)


class SimpleTARDISWorkflow(WorkflowLogging):
    show_progress_bars = Environment.allows_widget_display()
    enable_virtual_packet_logging = False
    log_level = None
    specific_log_level = None

    def __init__(self, configuration: Configuration, csvy: bool = False):
        """A simple TARDIS workflow that runs a simulation to convergence

        Parameters
        ----------
        configuration
            Configuration object for the simulation.
        csvy
            Set true if the configuration uses CSVY.
        """
        super().__init__(configuration, self.log_level, self.specific_log_level)
        atom_data = parse_atom_data(configuration)

        # set up states and solvers
        if csvy:
            self.simulation_state = SimulationState.from_csvy(configuration)
            assert np.isclose(
                self.simulation_state.v_inner_boundary.to(u.km / u.s).value,
                self.simulation_state.geometry.v_inner[0].to(u.km / u.s).value,
            ), (
                "If using csvy density input in the workflow, the initial v_inner_boundary must start at the first shell, see issue #3129."
            )

        else:
            self.simulation_state = SimulationState.from_config(
                configuration,
                atom_data=atom_data,
            )

        plasma_solver_factory = PlasmaSolverFactory(
            atom_data,
            configuration,
        )

        plasma_solver_factory.prepare_factory(
            self.simulation_state.abundance.index,
            "tardis.plasma.properties.property_collections",
            configuration,
        )

        self.plasma_solver = plasma_solver_factory.assemble(
            self.simulation_state.calculate_elemental_number_density(
                atom_data.atom_data.mass
            ),
            self.simulation_state.radiation_field_state,
            self.simulation_state.time_explosion,
            self.simulation_state._electron_densities,
        )

        line_interaction_type = configuration.plasma.line_interaction_type
        continuum_interactions = configuration.plasma.continuum_interaction

        self.opacity_solver = OpacitySolver(
            line_interaction_type,
            configuration.plasma.disable_line_scattering,
        )

        if line_interaction_type == "scatter":
            self.macro_atom_solver = None
        elif continuum_interactions.species:
            self.macro_atom_solver = ContinuumMacroAtomSolver(
                atom_data.levels,
                atom_data.lines,
                atom_data.photoionization_data,
                line_interaction_type=line_interaction_type,
            )
        else:
            self.macro_atom_solver = BoundBoundMacroAtomSolver(
                atom_data.levels,
                atom_data.lines,
                line_interaction_type,
            )
        self.transport_state = None
        self.transport_solver = MCTransportSolverClassic.from_config(
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

    def get_convergence_estimates(self) -> tuple[dict, object]:
        """Compute convergence estimates from the transport state

        Returns
        -------
        convergence_estimates
            Convergence estimates dictionary.
        estimated_radfield_properties
            Dilute radiation file and j_blues dataclass.
        """
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                self.transport_state.estimators_bulk,
                self.transport_state.estimators_line,
                self.transport_state.time_explosion,
                self.transport_state.time_of_simulation,
                self.transport_state.geometry_state.volume,
                self.transport_state.opacity_state.line_list_nu,
            )
        )

        estimated_t_radiative = estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        estimated_dilution_factor = estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor

        emitted_luminosity = calculate_filtered_luminosity(
            self.transport_state.emitted_packet_nu,
            self.transport_state.emitted_packet_luminosity,
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
        estimated_values: dict,
    ) -> bool:
        """Check convergence status for a dict of estimated values

        Parameters
        ----------
        estimated_values
            Estimates to check convergence.

        Returns
        -------
        converged
            If convergence has occurred.
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
        estimated_values: dict,
    ) -> dict:
        """Update the simulation state with new inputs computed from previous
        iteration estimates.

        Parameters
        ----------
        estimated_values
            Estimated from the previous iterations.

        Returns
        -------
        next_values
            Updated values for the simulation state.
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

    def solve_plasma(self, estimated_radfield_properties) -> None:
        """Update the plasma solution with the new radiation field estimates

        Parameters
        ----------
        estimated_radfield_properties
            The radiation field properties to use for updating the plasma.

        Raises
        ------
        ValueError
            If the plasma solver radiative rates type is unknown.
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
        if isinstance(self.macro_atom_solver, ContinuumMacroAtomSolver):
            continuum_property_solver = MCContinuumPropertiesSolver(
                self.plasma_solver.atomic_data
            )
            estimated_continuum_properties = continuum_property_solver.solve(
                self.transport_state.estimators_continuum,
                self.transport_state.time_of_simulation,
                self.transport_state.geometry_state.volume,
            )
            update_properties.update(
                gamma=estimated_continuum_properties.photo_ionization_rate_coefficient,
                alpha_stim_factor=estimated_continuum_properties.stimulated_recombination_rate_factor,
                bf_heating_coeff_estimator=self.transport_state.estimators_continuum.bf_heating_estimator,
                stim_recomb_cooling_coeff_estimator=self.transport_state.estimators_continuum.stim_recomb_cooling_estimator,
            )
        self.plasma_solver.update(**update_properties)

    def solve_opacity(self):
        """Solves the opacity state and any associated objects

        Returns
        -------
        dict
            opacity_state : tardis.opacities.opacity_state.OpacityState
                State of the line opacities
            macro_atom_state : tardis.opacities.macro_atom.macro_atom_state.MacroAtomState or None
                State of the macro atom
        """
        opacity_state = self.opacity_solver.solve(self.plasma_solver)

        if self.macro_atom_solver is None:
            macro_atom_state = None
        elif isinstance(self.macro_atom_solver, ContinuumMacroAtomSolver):
            macro_atom_state = self.macro_atom_solver.solve(
                self.plasma_solver.j_blues,
                opacity_state.beta_sobolev,
                self.plasma_solver.stimulated_emission_factor,
                self.plasma_solver.gamma_corr,
                self.plasma_solver.alpha_sp,
                self.plasma_solver.coll_deexc_coeff,
                self.plasma_solver.coll_exc_coeff,
                self.plasma_solver.electron_densities,
                self.plasma_solver.level_number_density,
                self.plasma_solver.delta_E_yg,
            )
            opacity_state.continuum_state.k_packet_idx = macro_atom_state.references_index.iloc[
                -1
            ]  # Hacky way to point to k-packet activation level - continuum state needs to be reexamined
        else:
            macro_atom_state = self.macro_atom_solver.solve(
                self.plasma_solver.j_blues,
                opacity_state.beta_sobolev,
                self.plasma_solver.stimulated_emission_factor,
            )

        return {
            "opacity_state": opacity_state,
            "macro_atom_state": macro_atom_state,
        }

    def solve_montecarlo(
        self,
        opacity_states: dict,
        no_of_real_packets: int,
        no_of_virtual_packets: int = 0,
    ) -> np.ndarray:
        """Solve the MonteCarlo process

        Parameters
        ----------
        opacity_states
            Opacity and (optionally) Macro Atom states.
        no_of_real_packets
            Number of real packets to simulate.
        no_of_virtual_packets
            Number of virtual packets to simulate per interaction.

        Returns
        -------
        virtual_packet_energies
            Array of unnormalized virtual packet energies in each frequency bin.
        """
        opacity_state = opacity_states["opacity_state"]
        macro_atom_state = opacity_states["macro_atom_state"]

        self.transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport_solver.run(
            self.transport_state,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = self.transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        return virtual_packet_energies

    def initialize_spectrum_solver(
        self,
        opacity_states: dict,
        virtual_packet_energies: np.ndarray | None = None,
    ) -> None:
        """Set up the spectrum solver

        Parameters
        ----------
        opacity_states
            Opacity and macro atom states.
        virtual_packet_energies
            Array of virtual packet energies binned by frequency.
        """
        # Set up spectrum solver
        self.spectrum_solver.transport_state = self.transport_state

        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_energies
            )

        if self.integrated_spectrum_settings is not None:
            # Set up spectrum solver integrator
            self.spectrum_solver.integrator_settings = (
                self.integrated_spectrum_settings
            )
            integrator_settings = self.spectrum_solver.integrator_settings
            formal_integrator = FormalIntegralSolver(
                integrator_settings.points,
                integrator_settings.interpolate_shells,
                getattr(integrator_settings, "method", None),
            )
            self.spectrum_solver.setup_optional_spectra(
                self.transport_state,
                virtual_packet_luminosity=None,
                integrator=formal_integrator,
                simulation_state=self.simulation_state,
                transport=self.transport_solver,
                plasma=self.plasma_solver,
                opacity_state=opacity_states["opacity_state"],
                macro_atom_state=opacity_states["macro_atom_state"],
            )

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        # Initialize iterations progress bar if showing progress bars
        if self.show_progress_bars:
            initialize_iterations_pbar(self.total_iterations)

        self.converged = False
        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                f"\n\tStarting iteration {(self.completed_iterations + 1):d} of {self.total_iterations:d}"
            )

            self.opacity_states = self.solve_opacity()

            virtual_packet_energies = self.solve_montecarlo(
                self.opacity_states, self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates()

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(estimated_radfield_properties)

            self.converged = self.check_convergence(estimated_values)
            self.completed_iterations += 1

            if self.converged and self.convergence_strategy.stop_if_converged:
                break

        if self.converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )
        self.opacity_states = self.solve_opacity()
        virtual_packet_energies = self.solve_montecarlo(
            self.opacity_states,
            self.final_iteration_packet_count,
            self.virtual_packet_count,
        )

        self.initialize_spectrum_solver(
            self.opacity_states,
            virtual_packet_energies,
        )
