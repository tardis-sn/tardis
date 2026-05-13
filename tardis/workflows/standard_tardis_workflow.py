import logging

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.simulation.base import PlasmaStateStorerMixin
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.util.environment import Environment
from tardis.visualization import ConvergencePlots
from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow

# logging support
logger = logging.getLogger(__name__)


class StandardTARDISWorkflow(
    SimpleTARDISWorkflow, PlasmaStateStorerMixin, HDFWriterMixin
):
    convergence_plots = None
    export_convergence_plots = False

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
        show_progress_bars=True,
        show_convergence_plots=False,
        convergence_plots_kwargs=None,
        csvy=False,
    ):
        if convergence_plots_kwargs is None:
            convergence_plots_kwargs = {}
        self.show_progress_bars = show_progress_bars
        self.log_level = log_level
        self.specific_log_level = specific_log_level
        self.enable_virtual_packet_logging = enable_virtual_packet_logging
        self.convergence_plots_kwargs = convergence_plots_kwargs

        SimpleTARDISWorkflow.__init__(self, configuration, csvy)

        # set up plasma storage
        PlasmaStateStorerMixin.__init__(
            self,
            iterations=self.total_iterations,
            no_of_shells=self.simulation_state.no_of_shells,
        )

        # Convergence plots
        if show_convergence_plots:
            (
                self.convergence_plots,
                self.export_convergence_plots,
            ) = self.initialize_convergence_plots()

    def initialize_convergence_plots(self):
        """Initialize the convergence plot attributes

        Returns
        -------
        ConvergencePlots
            The convergence plot instance
        bool
            If convergence plots are to be exported

        Raises
        ------
        RuntimeError
            Raised if run outside a notebook
        TypeError
            Raised if export_convergence_plots is not a bool
        """
        if not Environment.allows_widget_display():
            raise RuntimeError(
                "Convergence Plots cannot be displayed in command-line. Set show_convergence_plots "
                "to False."
            )

        convergence_plots = ConvergencePlots(
            iterations=self.total_iterations, **self.convergence_plots_kwargs
        )

        if "export_convergence_plots" in self.convergence_plots_kwargs:
            if not isinstance(
                self.convergence_plots_kwargs["export_convergence_plots"],
                bool,
            ):
                raise TypeError(
                    "Expected bool in export_convergence_plots argument"
                )
            export_convergence_plots = self.convergence_plots_kwargs[
                "export_convergence_plots"
            ]
        else:
            export_convergence_plots = False

        return convergence_plots, export_convergence_plots

    def get_convergence_estimates(self):
        """Compute convergence estimates from the transport state

        Returns
        -------
        dict
            Convergence estimates
        EstimatedRadiationFieldProperties
            Dilute radiation file and j_blues dataclass
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

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        self.converged = False
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

            self.opacity_states = self.solve_opacity()

            virtual_packet_energies = self.solve_montecarlo(
                self.opacity_states, self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates()

            if self.convergence_plots is not None:
                self.convergence_plots.update()

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
        self.store_plasma_state(
            self.completed_iterations,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_radiative,
            self.plasma_solver.electron_densities,
            self.simulation_state.t_inner,
        )
        self.reshape_plasma_state_store(self.completed_iterations)
        if self.convergence_plots is not None:
            self.get_convergence_estimates()
            self.convergence_plots.update(
                export_convergence_plots=self.export_convergence_plots,
                last=True,
            )
        self.initialize_spectrum_solver(
            self.opacity_states,
            virtual_packet_energies,
        )
