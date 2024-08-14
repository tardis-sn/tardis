import logging

from tardis.io.util import HDFWriterMixin
from tardis.simulation.base import PlasmaStateStorerMixin
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.util.base import is_notebook
from tardis.visualization import ConvergencePlots
from tardis.workflows.simple_simulation import SimpleSimulation

# logging support
logger = logging.getLogger(__name__)


class StandardSimulation(
    SimpleSimulation, PlasmaStateStorerMixin, HDFWriterMixin
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
        self.show_progress_bars = show_progress_bars
        self.log_level = log_level
        self.specific_log_level = specific_log_level
        self.enable_virtual_packet_logging = enable_virtual_packet_logging

        SimpleSimulation.__init__(self, configuration)

        # set up plasma storage
        PlasmaStateStorerMixin.__init__(
            self,
            iterations=self.total_iterations,
            no_of_shells=self.simulation_state.no_of_shells,
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

        estimated_t_radiative = estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        estimated_dilution_factor = estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor

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

        logger.info(
            f"\n\tLuminosity emitted   = {emitted_luminosity:.3e}\n"
            f"\tLuminosity absorbed  = {absorbed_luminosity:.3e}\n"
            f"\tLuminosity requested = {self.luminosity_requested:.3e}\n"
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

        if converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )
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
