import logging

from tardis.io.util import HDFWriterMixin
from tardis.simulation.base import PlasmaStateStorerMixin
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.util.base import is_notebook
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
        log_level="WARNING",
        specific_log_level=None,
        show_progress_bars=False,
        show_convergence_plots=False,
        convergence_plots_kwargs={},
    ):
        self.show_progress_bars = show_progress_bars
        self.log_level = log_level
        self.specific_log_level = specific_log_level
        self.enable_virtual_packet_logging = enable_virtual_packet_logging
        self.convergence_plots_kwargs = None
        self.show_convergence_plots = False

        SimpleTARDISWorkflow.__init__(self, configuration)

        # set up plasma storage
        PlasmaStateStorerMixin.__init__(
            self,
            iterations=self.total_iterations,
            no_of_shells=self.simulation_state.no_of_shells,
        )

        # Convergence plots
        show_convergence_plots = False
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
        if not is_notebook():
            raise RuntimeError(
                "Convergence Plots cannot be displayed in command-line. Set show_convergence_plots "
                "to False."
            )

        print("\n=== initialize_convergence_plots ===")
        print(f"total_iterations: {self.total_iterations}")
        print(f"convergence_plots_kwargs: {self.convergence_plots_kwargs}")

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
        print("\n=== get_convergence_estimates ===")
        print("Transport State:")
        # print(f"radfield_mc_estimators (first 10): {transport_state.radfield_mc_estimators[:10]}")
        # print(f"radfield_mc_estimators dtype: {transport_state.radfield_mc_estimators.dtype}")
        print(f"time_explosion: {transport_state.time_explosion}")
        print(f"time_of_simulation: {transport_state.time_of_simulation}")
        print(f"volume shape: {transport_state.geometry_state.volume.shape}")
        print(f"volume (first 10): {transport_state.geometry_state.volume[:10]}")

        print("\nSolving for estimated radfield properties:")
        print(f"line_list_nu (first 10): {transport_state.opacity_state.line_list_nu[:10]}")
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                transport_state.radfield_mc_estimators,
                transport_state.time_explosion,
                transport_state.time_of_simulation,
                transport_state.geometry_state.volume,
                transport_state.opacity_state.line_list_nu,
            )
        )
        print("\nRadfield Properties Result:")
        print(f"temperature (first 10): {estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature[:10]}")
        print(f"temperature shape: {estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature.shape}")

        estimated_t_radiative = estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        estimated_dilution_factor = estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor

        print("\nCalculating luminosities:")
        print(f"emitted_packet_nu (first 10): {transport_state.emitted_packet_nu[:10]}")
        print(f"emitted_packet_luminosity (first 10): {transport_state.emitted_packet_luminosity[:10]}")
        print(f"luminosity_nu_start: {self.luminosity_nu_start}")
        print(f"luminosity_nu_end: {self.luminosity_nu_end}")
        
        emitted_luminosity = calculate_filtered_luminosity(
            transport_state.emitted_packet_nu,
            transport_state.emitted_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )
        
        print(f"\nreabsorbed_packet_nu (first 10): {transport_state.reabsorbed_packet_nu[:10]}")
        print(f"reabsorbed_packet_luminosity (first 10): {transport_state.reabsorbed_packet_luminosity[:10]}")
        
        absorbed_luminosity = calculate_filtered_luminosity(
            transport_state.reabsorbed_packet_nu,
            transport_state.reabsorbed_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )

        print("\nCalculating luminosity ratios:")
        print(f"emitted_luminosity: {emitted_luminosity}")
        print(f"luminosity_requested: {self.luminosity_requested}")
        
        luminosity_ratios = (
            (emitted_luminosity / self.luminosity_requested).to(1).value
        )
        print(f"luminosity_ratios: {luminosity_ratios}")
        print(f"t_inner_update_exponent: {self.convergence_strategy.t_inner_update_exponent}")

        estimated_t_inner = (
            self.simulation_state.t_inner
            * luminosity_ratios
            ** self.convergence_strategy.t_inner_update_exponent
        )
        print(f"\nFinal estimated_t_inner: {estimated_t_inner}")

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

        print("\nEstimated Values:")
        print(f"t_radiative shape: {estimated_t_radiative.shape}")
        print(f"t_radiative (first 10): {estimated_t_radiative[:10]}")
        print(f"dilution_factor (first 10): {estimated_dilution_factor[:10]}")
        print(f"emitted_luminosity: {emitted_luminosity}")
        print(f"absorbed_luminosity: {absorbed_luminosity}")

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
        print("\n=== update_convergence_plot_data ===")
        print(f"Plot data keys: {plot_data_dict.keys()}")
        for name, (value, item_type) in plot_data_dict.items():
            print(f"\n{name}:")
            if hasattr(value, '__len__'):
                print(f"First 10 values: {value[:10]}")
                print(f"Shape: {getattr(value, 'shape', 'N/A')}")
                print(f"dtype: {getattr(value, 'dtype', 'N/A')}")
            else:
                print(f"Value: {value}")
            print(f"Type: {item_type}")
            self.convergence_plots.fetch_data(
                name=name,
                value=value,
                item_type=item_type,
            )

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        print("\n=== Starting run() ===")
        print(f"total_iterations: {self.total_iterations}")
        print(f"completed_iterations: {self.completed_iterations}")

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
            print(f"\n=== Iteration {self.completed_iterations + 1} ===")

            print("\nPlasma State:")
            print(f"dilution_factor (first 10): {self.simulation_state.dilution_factor[:10]}")
            print(f"t_radiative (first 10): {self.simulation_state.t_radiative[:10]}")
            print(f"electron_densities (first 10): {self.plasma_solver.electron_densities[:10]}")
            print(f"t_inner: {self.simulation_state.t_inner}")

            opacity_states = self.solve_opacity()
            print("\nOpacity States:")
            print(f"opacity_states type: {type(opacity_states)}")
            print(f"opacity_states keys: {getattr(opacity_states, '__dict__', {}).keys()}")

            transport_state, virtual_packet_energies = self.solve_montecarlo(
                opacity_states, self.real_packet_count
            )
            print("\nMonte Carlo Results:")
            print(f"transport_state type: {type(transport_state)}")
            print(f"virtual_packet_energies (first 10): {virtual_packet_energies[:10]}")

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
            opacity_states,
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
            self.get_convergence_estimates(transport_state)
            self.convergence_plots.update(
                export_convergence_plots=self.export_convergence_plots,
                last=True,
            )
        self.initialize_spectrum_solver(
            transport_state,
            opacity_states,
            virtual_packet_energies,
        )
