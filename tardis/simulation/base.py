import logging
import time
from collections import OrderedDict

import numpy as np
import pandas as pd
from astropy import units as u
from IPython.display import display

import tardis
from tardis import constants as const
from tardis.io.configuration.config_reader import ConfigurationError
from tardis.io.model.parse_atom_data import parse_atom_data
from tardis.io.model.parse_simulation_state import (
    parse_simulation_state,
)
from tardis.io.util import HDFWriterMixin
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.plasma.assembly.legacy_assembly import assemble_plasma
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral import FormalIntegrator
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
from tardis.transport.montecarlo.base import MonteCarloTransportSolver
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    MCContinuumPropertiesSolver,
)
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.macro_atom.macroatom_solver import MacroAtomSolver
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.util.base import is_notebook
from tardis.visualization import ConvergencePlots

logger = logging.getLogger(__name__)


class PlasmaStateStorerMixin:
    """Mixin class to provide the capability to the simulation object of
    storing plasma information and the inner boundary temperature during each
    MC iteration.

    Currently, storage for the dilution factor, the radiation temperature and
    the electron density in each cell is provided. Additionally, the
    temperature at the inner boundary is saved.
    """

    def __init__(self, iterations, no_of_shells):
        self.iterations_w = np.zeros((iterations, no_of_shells))
        self.iterations_t_rad = np.zeros((iterations, no_of_shells)) * u.K
        self.iterations_electron_densities = np.zeros(
            (iterations, no_of_shells)
        )
        self.iterations_t_inner = np.zeros(iterations) * u.K

    def store_plasma_state(self, i, w, t_rad, electron_densities, t_inner):
        """Store current plasma information and inner boundary temperature
        used in iterated i.

        Parameters
        ----------
        i : int
            current iteration index (0 for the first)
        w : np.ndarray
            dilution factor
        t_rad : astropy.units.Quantity
            radiation temperature
        electron_densities : np.ndarray
            electron density
        t_inner : astropy.units.Quantity
            temperature of inner boundary
        """
        self.iterations_w[i, :] = w
        self.iterations_t_rad[i, :] = t_rad
        self.iterations_electron_densities[i, :] = electron_densities.values
        self.iterations_t_inner[i] = t_inner

    def reshape_plasma_state_store(self, executed_iterations):
        """Reshapes the storage arrays in case convergence was reached before
        all specified iterations were executed.

        Parameters
        ----------
        executed_iterations : int
            iteration index, i.e. number of iterations executed minus one!
        """
        self.iterations_w = self.iterations_w[: executed_iterations + 1, :]
        self.iterations_t_rad = self.iterations_t_rad[
            : executed_iterations + 1, :
        ]
        self.iterations_electron_densities = self.iterations_electron_densities[
            : executed_iterations + 1, :
        ]
        self.iterations_t_inner = self.iterations_t_inner[
            : executed_iterations + 1
        ]


class Simulation(PlasmaStateStorerMixin, HDFWriterMixin):
    """A composite object containing all the required information for a
    simulation.

    Parameters
    ----------
    converged : bool
    iterations : int
    model : tardis.model.SimulationState
    plasma : tardis.plasma.BasePlasma
    transport : tardis.transport.montecarlo.MontecarloTransport
    opacity : tardis.opacities.opacity_solver.OpacitySolver
    macro_atom : tardis.opacities.macro_atom.macroatom_solver.MacroAtomSolver
    no_of_packets : int
    last_no_of_packets : int
    no_of_virtual_packets : int
    luminosity_nu_start : astropy.units.Quantity
    luminosity_nu_end : astropy.units.Quantity
    luminosity_requested : astropy.units.Quantity
    convergence_plots_kwargs: dict
    """

    hdf_properties = [
        "simulation_state",
        "plasma",
        "transport",
        "iterations_w",
        "iterations_t_rad",
        "iterations_electron_densities",
        "iterations_t_inner",
        "spectrum_solver",
    ]
    hdf_name = "simulation"

    def __init__(
        self,
        iterations,
        simulation_state,
        plasma,
        transport,
        opacity,
        macro_atom,
        no_of_packets,
        no_of_virtual_packets,
        luminosity_nu_start,
        luminosity_nu_end,
        last_no_of_packets,
        luminosity_requested,
        convergence_strategy,
        show_convergence_plots,
        convergence_plots_kwargs,
        show_progress_bars,
        spectrum_solver,
    ):
        super(Simulation, self).__init__(
            iterations, simulation_state.no_of_shells
        )

        self.converged = False
        self.iterations = iterations
        self.iterations_executed = 0
        self.simulation_state = simulation_state
        self.plasma = plasma
        self.transport = transport
        self.opacity = opacity
        self.macro_atom = macro_atom
        self.no_of_packets = no_of_packets
        self.last_no_of_packets = last_no_of_packets
        self.no_of_virtual_packets = no_of_virtual_packets
        self.luminosity_nu_start = luminosity_nu_start
        self.luminosity_nu_end = luminosity_nu_end
        self.luminosity_requested = luminosity_requested
        self.spectrum_solver = spectrum_solver
        self.show_progress_bars = show_progress_bars
        self.version = tardis.__version__

        # Convergence
        self.convergence_strategy = convergence_strategy
        self.converged = False
        self.consecutive_converges_count = 0

        # Convergence solvers
        self.t_rad_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_rad
        )
        self.w_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.w
        )
        self.t_inner_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_inner
        )

        if show_convergence_plots:
            if not is_notebook():
                raise RuntimeError(
                    "Convergence Plots cannot be displayed in command-line. Set show_convergence_plots "
                    "to False."
                )
            else:
                self.convergence_plots = ConvergencePlots(
                    iterations=self.iterations, **convergence_plots_kwargs
                )

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

        self._callbacks = OrderedDict()
        self._cb_next_id = 0

        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = (
            not self.plasma.continuum_interaction_species.empty
        )

    def estimate_t_inner(
        self,
        input_t_inner,
        luminosity_requested,
        emitted_luminosity,
        t_inner_update_exponent=-0.5,
    ):
        luminosity_ratios = (
            (emitted_luminosity / luminosity_requested).to(1).value
        )

        return input_t_inner * luminosity_ratios**t_inner_update_exponent

    def _get_convergence_status(
        self, t_rad, w, t_inner, estimated_t_rad, estimated_w, estimated_t_inner
    ):
        t_rad_converged = self.t_rad_convergence_solver.get_convergence_status(
            t_rad.value,
            estimated_t_rad.value,
            self.simulation_state.no_of_shells,
        )

        w_converged = self.w_convergence_solver.get_convergence_status(
            w, estimated_w, self.simulation_state.no_of_shells
        )

        t_inner_converged = (
            self.t_inner_convergence_solver.get_convergence_status(
                t_inner.value,
                estimated_t_inner.value,
                1,
            )
        )

        if np.all([t_rad_converged, w_converged, t_inner_converged]):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            # If an iteration has converged, require hold_iterations more
            # iterations to converge before we conclude that the Simulation
            # is converged.
            return self.consecutive_converges_count == hold_iterations + 1
        else:
            self.consecutive_converges_count = 0
            return False

    def advance_state(self, emitted_luminosity):
        """
        Advances the state of the model and the plasma for the next
        iteration of the simulation. Returns True if the convergence criteria
        are met, else False.

        Returns
        -------
            converged : bool
        """
        estimated_radfield_properties = (
            self.transport.radfield_prop_solver.solve(
                self.transport.transport_state.radfield_mc_estimators,
                self.transport.transport_state.time_explosion,
                self.transport.transport_state.time_of_simulation,
                self.transport.transport_state.geometry_state.volume,
                self.transport.transport_state.opacity_state.line_list_nu,
            )
        )

        estimated_t_rad = (
            estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        )
        estimated_dilution_factor = (
            estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor
        )

        estimated_t_inner = self.estimate_t_inner(
            self.simulation_state.t_inner,
            self.luminosity_requested,
            emitted_luminosity,
            t_inner_update_exponent=self.convergence_strategy.t_inner_update_exponent,
        )

        converged = self._get_convergence_status(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_rad,
            estimated_dilution_factor,
            estimated_t_inner,
        )

        # calculate_next_plasma_state equivalent
        next_t_radiative = self.t_rad_convergence_solver.converge(
            self.simulation_state.t_radiative,
            estimated_t_rad,
        )
        next_dilution_factor = self.w_convergence_solver.converge(
            self.simulation_state.dilution_factor,
            estimated_dilution_factor,
        )
        if (
            self.iterations_executed + 1
        ) % self.convergence_strategy.lock_t_inner_cycles == 0:
            next_t_inner = self.t_inner_convergence_solver.converge(
                self.simulation_state.t_inner,
                estimated_t_inner,
            )
        else:
            next_t_inner = self.simulation_state.t_inner

        if hasattr(self, "convergence_plots"):
            self.convergence_plots.fetch_data(
                name="t_inner",
                value=self.simulation_state.t_inner.value,
                item_type="value",
            )
            self.convergence_plots.fetch_data(
                name="t_rad",
                value=self.simulation_state.t_radiative,
                item_type="iterable",
            )
            self.convergence_plots.fetch_data(
                name="w",
                value=self.simulation_state.dilution_factor,
                item_type="iterable",
            )
            self.convergence_plots.fetch_data(
                name="velocity",
                value=self.simulation_state.velocity,
                item_type="iterable",
            )

        self.log_plasma_state(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            next_t_radiative,
            next_dilution_factor,
            next_t_inner,
        )
        self.simulation_state.t_radiative = next_t_radiative
        self.simulation_state.dilution_factor = next_dilution_factor
        self.simulation_state.blackbody_packet_source.temperature = next_t_inner

        radiation_field = DilutePlanckianRadiationField(
            temperature=self.simulation_state.t_radiative,
            dilution_factor=self.simulation_state.dilution_factor,
        )
        update_properties = dict(
            dilute_planckian_radiation_field=radiation_field
        )

        # model.calculate_j_blues() equivalent
        # model.update_plasmas() equivalent
        # Bad test to see if this is a nlte run

        if "nlte_data" in self.plasma.outputs_dict:
            self.plasma.store_previous_properties()

        # JBlues solver
        if (
            self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "blackbody"
        ):
            planckian_radiation_field = (
                radiation_field.to_planckian_radiation_field()
            )
            j_blues = planckian_radiation_field.calculate_mean_intensity(
                self.plasma.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma.atomic_data.lines.index
            )
        elif (
            self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "dilute-blackbody"
        ):
            j_blues = radiation_field.calculate_mean_intensity(
                self.plasma.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma.atomic_data.lines.index
            )
        elif (
            self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "detailed"
        ):
            update_properties["j_blues"] = pd.DataFrame(
                estimated_radfield_properties.j_blues,
                index=self.plasma.atomic_data.lines.index,
            )
        else:
            raise ValueError(
                f"radiative_rates_type type unknown - {self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE}"
            )

        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.

        radfield_mc_estimators = (
            self.transport.transport_state.radfield_mc_estimators
        )

        if "gamma" in self.plasma.outputs_dict:
            continuum_property_solver = MCContinuumPropertiesSolver(
                self.atom_data
            )
            estimated_continuum_properties = continuum_property_solver.solve(
                radfield_mc_estimators,
                self.transport.transport_state.time_of_simulation,
                self.transport.transport_state.geometry_state.volume,
            )
            update_properties.update(
                gamma=estimated_continuum_properties.photo_ion_coeff,
                alpha_stim_coeff=estimated_continuum_properties.stim_recomb_estimator,
                bf_heating_coeff_estimator=radfield_mc_estimators.bf_heating_estimator,
                stim_recomb_cooling_coeff_estimator=radfield_mc_estimators.stim_recomb_cooling_estimator,
            )

        self.plasma.update(**update_properties)

        return converged

    def iterate(self, no_of_packets, no_of_virtual_packets=0):
        logger.info(
            f"\n\tStarting iteration {(self.iterations_executed + 1):d} of {self.iterations:d}"
        )

        opacity_state = self.opacity.solve(self.plasma)
        if self.macro_atom is not None:
            if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
                macro_atom_state = MacroAtomState.from_legacy_plasma(
                    self.plasma
                )  # TODO: Impliment
            else:
                macro_atom_state = self.macro_atom.solve(
                    self.plasma,
                    self.plasma.atomic_data,
                    opacity_state.tau_sobolev,
                    self.plasma.stimulated_emission_factor,
                )

        transport_state = self.transport.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            self.plasma,
            no_of_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.iterations_executed,
        )

        v_packets_energy_hist = self.transport.run(
            transport_state,
            iteration=self.iterations_executed,
            total_iterations=self.iterations,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = (
            self.transport.transport_state.packet_collection.output_energies
        )
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        emitted_luminosity = calculate_filtered_luminosity(
            transport_state.emitted_packet_nu,
            transport_state.emitted_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )
        reabsorbed_luminosity = calculate_filtered_luminosity(
            transport_state.reabsorbed_packet_nu,
            transport_state.reabsorbed_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )
        if hasattr(self, "convergence_plots"):
            self.convergence_plots.fetch_data(
                name="Emitted",
                value=emitted_luminosity.value,
                item_type="value",
            )
            self.convergence_plots.fetch_data(
                name="Absorbed",
                value=reabsorbed_luminosity.value,
                item_type="value",
            )
            self.convergence_plots.fetch_data(
                name="Requested",
                value=self.luminosity_requested.value,
                item_type="value",
            )

        self.log_run_results(emitted_luminosity, reabsorbed_luminosity)
        self.iterations_executed += 1
        return emitted_luminosity, v_packets_energy_hist

    def run_convergence(self):
        """
        run the simulation
        """
        start_time = time.time()
        while self.iterations_executed < self.iterations - 1:
            self.store_plasma_state(
                self.iterations_executed,
                self.simulation_state.dilution_factor,
                self.simulation_state.t_radiative,
                self.plasma.electron_densities,
                self.simulation_state.t_inner,
            )
            emitted_luminosity, v_packets_energy_hist = self.iterate(
                self.no_of_packets
            )
            self.converged = self.advance_state(emitted_luminosity)
            if hasattr(self, "convergence_plots"):
                self.convergence_plots.update()
            self._call_back()
            if self.converged:
                if self.convergence_strategy.stop_if_converged:
                    break

        logger.info(
            f"\n\tSimulation finished in {self.iterations_executed:d} iterations "
            f"\n\tSimulation took {(time.time() - start_time):.2f} s\n"
        )

    def run_final(self):
        """
        run the last iteration of the simulation
        """
        self.store_plasma_state(
            self.iterations_executed,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_radiative,
            self.plasma.electron_densities,
            self.simulation_state.t_inner,
        )
        emitted_luminosity, v_packets_energy_hist = self.iterate(
            self.last_no_of_packets, self.no_of_virtual_packets
        )

        # Set up spectrum solver integrator and virtual spectrum
        self.spectrum_solver.setup_optional_spectra(
            self.transport.transport_state,
            v_packets_energy_hist,
            FormalIntegrator(
                self.simulation_state, self.plasma, self.transport
            ),
        )

        self.reshape_plasma_state_store(self.iterations_executed)
        if hasattr(self, "convergence_plots"):
            self.convergence_plots.fetch_data(
                name="t_inner",
                value=self.simulation_state.t_inner.value,
                item_type="value",
            )
            self.convergence_plots.update(
                export_convergence_plots=self.export_convergence_plots,
                last=True,
            )

        self._call_back()

    def log_plasma_state(
        self,
        t_rad,
        dilution_factor,
        t_inner,
        next_t_rad,
        next_dilution_factor,
        next_t_inner,
        log_sampling=5,
    ):
        """
        Logging the change of the plasma state

        Parameters
        ----------
        t_rad : astropy.units.Quanity
            current t_rad
        dilution_factor : np.ndarray
            current dilution_factor
        next_t_rad : astropy.units.Quanity
            next t_rad
        next_dilution_factor : np.ndarray
            next dilution_factor
        log_sampling : int
            the n-th shells to be plotted

        Returns
        -------
        """
        plasma_state_log = pd.DataFrame(
            index=np.arange(len(t_rad)),
            columns=["t_rad", "next_t_rad", "w", "next_w"],
        )
        plasma_state_log["t_rad"] = t_rad
        plasma_state_log["next_t_rad"] = next_t_rad
        plasma_state_log["w"] = dilution_factor
        plasma_state_log["next_w"] = next_dilution_factor
        plasma_state_log.columns.name = "Shell No."

        if is_notebook():
            logger.info("\n\tPlasma stratification:")

            # Displaying the DataFrame only when the logging level is NOTSET, DEBUG or INFO
            if logger.level <= logging.INFO:
                if not logger.filters:
                    display(
                        plasma_state_log.iloc[::log_sampling].style.format(
                            "{:.3g}"
                        )
                    )
                elif logger.filters[0].log_level == 20:
                    display(
                        plasma_state_log.iloc[::log_sampling].style.format(
                            "{:.3g}"
                        )
                    )
        else:
            output_df = ""
            plasma_output = plasma_state_log.iloc[::log_sampling].to_string(
                float_format=lambda x: f"{x:.3g}",
                justify="center",
            )
            for value in plasma_output.split("\n"):
                output_df = output_df + f"\t{value}\n"
            logger.info("\n\tPlasma stratification:")
            logger.info(f"\n{output_df}")

        logger.info(
            f"\n\tCurrent t_inner = {t_inner:.3f}\n\tExpected t_inner for next iteration = {next_t_inner:.3f}\n"
        )

    def log_run_results(self, emitted_luminosity, absorbed_luminosity):
        logger.info(
            f"\n\tLuminosity emitted   = {emitted_luminosity:.3e}\n"
            f"\tLuminosity absorbed  = {absorbed_luminosity:.3e}\n"
            f"\tLuminosity requested = {self.luminosity_requested:.3e}\n"
        )

    def _call_back(self):
        for cb, args in self._callbacks.values():
            cb(self, *args)

    def add_callback(self, cb_func, *args):
        """
        Add a function which will be called
        after every iteration.

        The cb_func signature must look like:
        cb_func(simulation, extra_arg1, ...)

        Parameters
        ----------
        cb_func : callable
            The callback function
        arg1 :
            The first additional arguments passed to the callable function
        ...

        Returns
        -------
        : int
            The callback ID
        """
        cb_id = self._cb_next_id
        self._callbacks[cb_id] = (cb_func, args)
        self._cb_next_id += 1
        return cb_id

    def remove_callback(self, id):
        """
        Remove the callback with a specific ID (which was returned by
        add_callback)

        Parameters
        ----------
        id : int
            The callback ID

        Returns
        -------
        : True if the callback was successfully removed.
        """
        try:
            del self._callbacks[id]
            return True
        except KeyError:
            logger.debug(f"Call Back was not found in {self._callbacks.keys()}")
            return False

    @classmethod
    def from_config(
        cls,
        config,
        packet_source=None,
        virtual_packet_logging=False,
        show_convergence_plots=False,
        show_progress_bars=True,
        legacy_mode_enabled=False,
        atom_data=None,
        plasma=None,
        transport=None,
        opacity=None,
        macro_atom=None,
        **kwargs,
    ):
        """
        Create a simulation instance from the provided configuration.

        Parameters
        ----------
        config : object
            The configuration object for the simulation.
        packet_source : object, optional
            The packet source for the simulation.
        virtual_packet_logging : bool, optional
            Flag indicating virtual packet logging.
        show_convergence_plots : bool, optional
            Flag indicating whether to show convergence plots.
        show_progress_bars : bool, optional
            Flag indicating whether to show progress bars.
        legacy_mode_enabled : bool, optional
            Flag indicating if legacy mode is enabled.
        atom_data : object, optional
            The atom data for the simulation.
        plasma : object, optional
            The plasma object for the simulation.
        transport : object, optional
            The transport solver for the simulation.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        object
            The created simulation instance.
        """
        # Allow overriding some config structures. This is useful in some
        # unit tests, and could be extended in all the from_config classmethods.

        atom_data = parse_atom_data(config, atom_data=atom_data)
        simulation_state = parse_simulation_state(
            config, packet_source, legacy_mode_enabled, kwargs, atom_data
        )
        if plasma is None:
            plasma = assemble_plasma(
                config,
                simulation_state,
                atom_data=atom_data,
            )

        if (transport is not None) and (packet_source is not None):
            raise ConfigurationError(
                "Cannot specify packet_source and transport at the same time."
            )
        if transport is None:
            transport = MonteCarloTransportSolver.from_config(
                config,
                packet_source=simulation_state.packet_source,
                enable_virtual_packet_logging=virtual_packet_logging,
            )
        if opacity is None:
            opacity = OpacitySolver(
                config.plasma.line_interaction_type,
                config.plasma.disable_line_scattering,
            )
        if macro_atom is None:
            if config.plasma.line_interaction_type in (
                "downbranch",
                "macroatom",
            ):
                macro_atom = MacroAtomSolver()

        convergence_plots_config_options = [
            "plasma_plot_config",
            "t_inner_luminosities_config",
            "plasma_cmap",
            "t_inner_luminosities_colors",
            "export_convergence_plots",
        ]
        convergence_plots_kwargs = {}
        for item in set(convergence_plots_config_options).intersection(
            kwargs.keys()
        ):
            convergence_plots_kwargs[item] = kwargs[item]

        luminosity_nu_start = config.supernova.luminosity_wavelength_end.to(
            u.Hz, u.spectral()
        )

        if u.isclose(
            config.supernova.luminosity_wavelength_start, 0 * u.angstrom
        ):
            luminosity_nu_end = np.inf * u.Hz
        else:
            luminosity_nu_end = (
                const.c / config.supernova.luminosity_wavelength_start
            ).to(u.Hz)

        last_no_of_packets = config.montecarlo.last_no_of_packets
        if last_no_of_packets is None or last_no_of_packets < 0:
            last_no_of_packets = config.montecarlo.no_of_packets
        last_no_of_packets = int(last_no_of_packets)

        spectrum_solver = SpectrumSolver.from_config(config)

        return cls(
            iterations=config.montecarlo.iterations,
            simulation_state=simulation_state,
            plasma=plasma,
            transport=transport,
            opacity=opacity,
            macro_atom=macro_atom,
            show_convergence_plots=show_convergence_plots,
            no_of_packets=int(config.montecarlo.no_of_packets),
            no_of_virtual_packets=int(config.montecarlo.no_of_virtual_packets),
            luminosity_nu_start=luminosity_nu_start,
            luminosity_nu_end=luminosity_nu_end,
            last_no_of_packets=last_no_of_packets,
            luminosity_requested=config.supernova.luminosity_requested.cgs,
            convergence_strategy=config.montecarlo.convergence_strategy,
            convergence_plots_kwargs=convergence_plots_kwargs,
            show_progress_bars=show_progress_bars,
            spectrum_solver=spectrum_solver,
        )
