import logging

from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.io.configuration.config_reader import Configuration
from tardis.model.geometry.radial1d_nonhomologous import (
    NonhomologousRadial1DGeometry,
)
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
)
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.modes.nonhomologous.opacity_solver import (
    OpacitySolver,
)
from tardis.transport.montecarlo.modes.nonhomologous.plasma_assembly_base import (
    PlasmaSolverFactory,
)
from tardis.transport.montecarlo.modes.nonhomologous.solver import (
    MCTransportSolverNonhomologous,
)
from tardis.util.environment import Environment
from tardis.visualization import ConvergencePlots
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow

# logging support
logger = logging.getLogger(__name__)


class NonhomologousTARDISWorkflow(StandardTARDISWorkflow):
    def __init__(
        self,
        configuration: Configuration,
        csvy: bool = False,
        log_level: str | None = None,
        specific_log_level: bool | None = None,
        show_convergence_plots: bool = False,
    ):
        """
        Inherits from StandardTARDISWorkflow and overrides the components that
        differ for non-homologous expansion: the geometry, plasma solver,
        opacity solver, and transport solver.

        Parameters
        ----------
        configuration
            Configuration object for the simulation.
        csvy
            Set true if the configuration uses CSVY.
        log_level
            Sets the logging Level for the logger
        specific_log_level
            Allows for logging on the specified logging levels
        show_convergence_plots
            Whether to display convergence plots while iterating
        """
        super().__init__(
            configuration,
            csvy=csvy,
            log_level=log_level,
            specific_log_level=specific_log_level,
            show_convergence_plots=show_convergence_plots,
        )
        atom_data = parse_atom_data(configuration)

        # Replace the default geometry of the SimpleTARDISWorkflow
        geometry = self.simulation_state.geometry
        t_exp = self.simulation_state.time_explosion
        self.simulation_state.geometry = NonhomologousRadial1DGeometry(
            r_inner=geometry.v_inner * t_exp,
            r_outer=geometry.v_outer * t_exp,
            v_inner=geometry.v_inner,
            v_outer=geometry.v_outer,
            r_inner_boundary=geometry.v_inner_boundary * t_exp,
            r_outer_boundary=geometry.v_outer_boundary * t_exp,
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
            self.simulation_state._electron_densities,
        )

        line_interaction_type = configuration.plasma.line_interaction_type
        self.opacity_solver = OpacitySolver(
            self.simulation_state.geometry.velocity_gradient,
            line_interaction_type,
            configuration.plasma.disable_line_scattering,
        )

        # TODO: continuum support
        if line_interaction_type == "scatter":
            self.macro_atom_solver = None
        else:
            self.macro_atom_solver = BoundBoundMacroAtomSolver(
                atom_data.levels,
                atom_data.lines,
                line_interaction_type,
            )

        self.transport_solver = MCTransportSolverNonhomologous.from_config(
            configuration,
            packet_source=self.simulation_state.packet_source,
            enable_virtual_packet_logging=self.enable_virtual_packet_logging,
        )

    def get_convergence_estimates(self) -> tuple[dict, object]:
        """Compute convergence estimates from the transport state

        Returns
        -------
        convergence_estimates
            Convergence estimates dictionary.
        estimated_radfield_properties
            Dilute radiation field and j_blues dataclass.
        """
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                self.transport_state.estimators_bulk,
                self.transport_state.estimators_line,
                self.transport_state.geometry_state.velocity_gradient,
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
        absorbed_luminosity = calculate_filtered_luminosity(
            self.transport_state.reabsorbed_packet_nu,
            self.transport_state.reabsorbed_packet_luminosity,
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

        logger.info(
            f"\n\tLuminosity emitted   = {emitted_luminosity:.3e}\n"
            f"\tLuminosity absorbed  = {absorbed_luminosity:.3e}\n"
            f"\tLuminosity requested = {self.luminosity_requested:.3e}\n"
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
