import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.model.geometry.radial1d_nonhomologous import (
    NonhomologousRadial1DGeometry,
)
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
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
from tardis.transport.montecarlo.modes.nonhomologous.solver import (
    MCTransportSolverNonhomologous,
)
from tardis.transport.montecarlo.progress_bars import initialize_iterations_pbar
from tardis.util.environment import Environment
from tardis.workflows.workflow_logging import WorkflowLogging

# logging support
logger = logging.getLogger(__name__)


class NonhomologousTARDISWorkflow(WorkflowLogging):
    show_progress_bars = Environment.allows_widget_display()
    enable_virtual_packet_logging = False
    log_level = None
    specific_log_level = None

    def __init__(self, configuration: Configuration, csvy: bool = False):
        super().__init__(configuration, self.log_level, self.specific_log_level)
        atom_data = parse_atom_data(configuration)

        # Set up the simulation state
        if csvy:
            self.simulation_state = SimulationState.from_csvy(configuration)
            assert np.isclose(
                self.simulation_state.v_inner_boundary.to(u.km / u.s).value,
                self.simulation_state.geometry.v_inner[0].to(u.km / u.s).value,
            ), (
                "If using csvy density input in the workflow, the initial "
                "v_inner_boundary must start at the first shell, see issue #3129."
            )
        else:
            self.simulation_state = SimulationState.from_config(
                configuration,
                atom_data=atom_data,
            )

        # Replace homologous geometry in the simulation state with nonhomologous geometry object
        geometry = self.simulation_state.geometry
        t_exp = self.simulation_state.time_explosion
        self.simulation_state.geometry = NonhomologousRadial1DGeometry(
            r_inner=geometry.v_inner * t_exp,
            r_outer=geometry.v_outer * t_exp,
            v_inner=geometry.v_inner,
            v_outer=geometry.v_outer,
            r_inner_boundary=None,
            r_outer_boundary=None,
            v_inner_boundary=geometry.v_inner_boundary,
            v_outer_boundary=geometry.v_outer_boundary,
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

        # Should always be macroatom for now
        line_interaction_type = configuration.plasma.line_interaction_type

        self.opacity_solver = OpacitySolver(
            line_interaction_type,
            configuration.plasma.disable_line_scattering,
        )

        self.macro_atom_solver = BoundBoundMacroAtomSolver(
            atom_data.levels,
            atom_data.lines,
            line_interaction_type,
        )

        # Solve the opacity state once at the beginning
        opacity_state = self.opacity_solver.solve(self.plasma_solver)
        macro_atom_state = self.macro_atom_solver.solve(
            self.plasma_solver.j_blues,
            opacity_state.beta_sobolev,
            self.plasma_solver.stimulated_emission_factor,
        )
        self.opacity_states = {
            "opacity_state": opacity_state,
            "macro_atom_state": macro_atom_state,
        }

        # Initialize the transport state and solver
        # Disable full relativity (not supported yet for nonhomology)
        self.transport_state = None
        self.transport_solver = MCTransportSolverNonhomologous.from_config(
            configuration,
            packet_source=self.simulation_state.packet_source,
            enable_virtual_packet_logging=self.enable_virtual_packet_logging,
        )

        self.transport_solver.montecarlo_configuration.ENABLE_FULL_RELATIVITY = (
            False
        )

        # connor-mcclellan: the rest of the init is unmodified from the simple workflow

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

        # Monte Carlo settings
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

        # Spectrum settings
        self.integrated_spectrum_settings = configuration.spectrum.integrated
        self.spectrum_solver = SpectrumSolver.from_config(configuration)


    def solve_montecarlo(
        self,
        no_of_real_packets: int,
        no_of_virtual_packets: int = 0,
    ) -> np.ndarray:
        """Run non-homologous Monte Carlo transport.

        Parameters
        ----------
        no_of_real_packets : int
            Number of real packets to simulate.
        no_of_virtual_packets : int, optional
            Number of virtual packets per interaction, by default 0.

        Returns
        -------
        virtual_packet_energies : np.ndarray
            Array of unnormalized virtual packet energies in each frequency bin.
        """
        opacity_state = self.opacity_states["opacity_state"]
        macro_atom_state = self.opacity_states["macro_atom_state"]

        self.transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=0,
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
        virtual_packet_energies: np.ndarray | None = None,
    ) -> None:
        """Set up the spectrum solver

        Parameters
        ----------
        virtual_packet_energies : np.ndarray, optional
            Array of virtual packet energies binned by frequency.
        """
        self.spectrum_solver.transport_state = self.transport_state

        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_energies
            )

        if self.integrated_spectrum_settings is not None:
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
                opacity_state=self.opacity_states["opacity_state"],
                macro_atom_state=self.opacity_states["macro_atom_state"],
            )

    def run(self) -> None:
        """Run the non-homologous TARDIS simulation - single iteration for now."""
        if self.show_progress_bars:
            initialize_iterations_pbar(self.total_iterations)

        virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count,
            self.virtual_packet_count,
        )

        self.initialize_spectrum_solver(virtual_packet_energies)
