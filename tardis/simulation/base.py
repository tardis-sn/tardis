import logging
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from IPython.display import display

import tardis
from tardis import constants as const
from tardis.io.atom_data.base import AtomData
from tardis.io.configuration.config_reader import ConfigurationError
from tardis.io.util import HDFWriterMixin
from tardis.model import SimulationState
from tardis.model.parse_input import initialize_packet_source
from tardis.montecarlo.base import MonteCarloTransportSolver
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.util.base import is_notebook
from tardis.visualization import ConvergencePlots

# Adding logging support
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
    transport : tardis.montecarlo.MontecarloTransport
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
    ]
    hdf_name = "simulation"

    def __init__(
        self,
        iterations,
        simulation_state,
        plasma,
        transport,
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
        self.no_of_packets = no_of_packets
        self.last_no_of_packets = last_no_of_packets
        self.no_of_virtual_packets = no_of_virtual_packets
        self.luminosity_nu_start = luminosity_nu_start
        self.luminosity_nu_end = luminosity_nu_end
        self.luminosity_requested = luminosity_requested
        self.show_progress_bars = show_progress_bars
        self.version = tardis.__version__

        if convergence_strategy.type in ("damped"):
            self.convergence_strategy = convergence_strategy
            self.converged = False
            self.consecutive_converges_count = 0
        elif convergence_strategy.type in ("custom"):
            raise NotImplementedError(
                "Convergence strategy type is custom; "
                "you need to implement your specific treatment!"
            )
        else:
            raise ValueError(
                f"Convergence strategy type is "
                f"not damped or custom "
                f"- input is {convergence_strategy.type}"
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

        self.transport.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED = (
            not self.plasma.continuum_interaction_species.empty
        )

    def estimate_t_inner(
        self, input_t_inner, luminosity_requested, t_inner_update_exponent=-0.5
    ):
        emitted_luminosity = (
            self.transport.transport_state.calculate_emitted_luminosity(
                self.luminosity_nu_start, self.luminosity_nu_end
            )
        )

        luminosity_ratios = (
            (emitted_luminosity / luminosity_requested).to(1).value
        )

        return input_t_inner * luminosity_ratios**t_inner_update_exponent

    @staticmethod
    def damped_converge(value, estimated_value, damping_factor):
        # FIXME: Should convergence strategy have its own class containing this
        # as a method
        return value + damping_factor * (estimated_value - value)

    def _get_convergence_status(
        self, t_rad, w, t_inner, estimated_t_rad, estimated_w, estimated_t_inner
    ):
        # FIXME: Move the convergence checking in its own class.
        no_of_shells = self.simulation_state.no_of_shells

        convergence_t_rad = (
            abs(t_rad - estimated_t_rad) / estimated_t_rad
        ).value
        convergence_w = abs(w - estimated_w) / estimated_w
        convergence_t_inner = (
            abs(t_inner - estimated_t_inner) / estimated_t_inner
        ).value

        fraction_t_rad_converged = (
            np.count_nonzero(
                convergence_t_rad < self.convergence_strategy.t_rad.threshold
            )
            / no_of_shells
        )

        t_rad_converged = (
            fraction_t_rad_converged > self.convergence_strategy.fraction
        )

        fraction_w_converged = (
            np.count_nonzero(
                convergence_w < self.convergence_strategy.w.threshold
            )
            / no_of_shells
        )

        w_converged = fraction_w_converged > self.convergence_strategy.fraction

        t_inner_converged = (
            convergence_t_inner < self.convergence_strategy.t_inner.threshold
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

    def advance_state(self):
        """
        Advances the state of the model and the plasma for the next
        iteration of the simulation. Returns True if the convergence criteria
        are met, else False.

        Returns
        -------
            converged : bool
        """
        (
            estimated_t_rad,
            estimated_dilution_factor,
        ) = self.transport.transport_state.calculate_radiationfield_properties()
        estimated_t_inner = self.estimate_t_inner(
            self.simulation_state.t_inner,
            self.luminosity_requested,
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
        # FIXME: Should convergence strategy have its own class?
        next_t_radiative = self.damped_converge(
            self.simulation_state.t_radiative,
            estimated_t_rad,
            self.convergence_strategy.t_rad.damping_constant,
        )
        next_dilution_factor = self.damped_converge(
            self.simulation_state.dilution_factor,
            estimated_dilution_factor,
            self.convergence_strategy.w.damping_constant,
        )
        if (
            self.iterations_executed + 1
        ) % self.convergence_strategy.lock_t_inner_cycles == 0:
            next_t_inner = self.damped_converge(
                self.simulation_state.t_inner,
                estimated_t_inner,
                self.convergence_strategy.t_inner.damping_constant,
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

        # model.calculate_j_blues() equivalent
        # model.update_plasmas() equivalent
        # Bad test to see if this is a nlte run
        if "nlte_data" in self.plasma.outputs_dict:
            self.plasma.store_previous_properties()

        update_properties = dict(
            t_rad=self.simulation_state.t_radiative,
            w=self.simulation_state.dilution_factor,
        )
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.

        estimators = self.transport.transport_state.radfield_mc_estimators
        if "j_blue_estimator" in self.plasma.outputs_dict:
            update_properties.update(
                t_inner=next_t_inner,
                j_blue_estimator=estimators.j_blue_estimator,
            )
        if "gamma_estimator" in self.plasma.outputs_dict:
            update_properties.update(
                gamma_estimator=estimators.photo_ion_estimator,
                alpha_stim_estimator=estimators.stim_recomb_estimator,
                bf_heating_coeff_estimator=estimators.bf_heating_estimator,
                stim_recomb_cooling_coeff_estimator=estimators.stim_recomb_cooling_estimator,
            )

        self.plasma.update(**update_properties)

        return converged

    def iterate(self, no_of_packets, no_of_virtual_packets=0):
        logger.info(
            f"\n\tStarting iteration {(self.iterations_executed + 1):d} of {self.iterations:d}"
        )

        transport_state = self.transport.initialize_transport_state(
            self.simulation_state,
            self.plasma,
            no_of_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.iterations_executed,
        )

        self.transport.run(
            transport_state,
            time_explosion=self.simulation_state.time_explosion,
            iteration=self.iterations_executed,
            total_iterations=self.iterations,
            show_progress_bars=self.show_progress_bars,
        )

        output_energy = (
            self.transport.transport_state.packet_collection.output_energies
        )
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        emitted_luminosity = (
            self.transport.transport_state.calculate_emitted_luminosity(
                self.luminosity_nu_start, self.luminosity_nu_end
            )
        )
        reabsorbed_luminosity = (
            self.transport.transport_state.calculate_reabsorbed_luminosity(
                self.luminosity_nu_start, self.luminosity_nu_end
            )
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
            self.iterate(self.no_of_packets)
            self.converged = self.advance_state()
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
        self.iterate(self.last_no_of_packets, self.no_of_virtual_packets)

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
        w : astropy.units.Quanity
            current w
        next_t_rad : astropy.units.Quanity
            next t_rad
        next_w : astropy.units.Quanity
            next_w
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
        **kwargs,
    ):
        """
        Create a new Simulation instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        **kwargs
            Allow overriding some structures, such as model, plasma, atomic data
            and the transport, instead of creating them from the configuration
            object.

        Returns
        -------
        Simulation
        """
        # Allow overriding some config structures. This is useful in some
        # unit tests, and could be extended in all the from_config classmethods.

        atom_data = kwargs.get("atom_data", None)
        if atom_data is None:
            if "atom_data" in config:
                if Path(config.atom_data).is_absolute():
                    atom_data_fname = Path(config.atom_data)
                else:
                    atom_data_fname = (
                        Path(config.config_dirname) / config.atom_data
                    )

            else:
                raise ValueError(
                    "No atom_data option found in the configuration."
                )

            logger.info(f"\n\tReading Atomic Data from {atom_data_fname}")

            try:
                atom_data = AtomData.from_hdf(atom_data_fname)
            except TypeError as e:
                print(
                    e,
                    "Error might be from the use of an old-format of the atomic database, \n"
                    "please see https://github.com/tardis-sn/tardis-refdata/tree/master/atom_data"
                    " for the most recent version.",
                )
                raise
        if "model" in kwargs:
            simulation_state = kwargs["model"]
        else:
            if hasattr(config, "csvy_model"):
                simulation_state = SimulationState.from_csvy(
                    config,
                    atom_data=atom_data,
                    legacy_mode_enabled=legacy_mode_enabled,
                )
            else:
                simulation_state = SimulationState.from_config(
                    config,
                    atom_data=atom_data,
                    legacy_mode_enabled=legacy_mode_enabled,
                )
            if packet_source is not None:
                simulation_state.packet_source = initialize_packet_source(
                    config,
                    simulation_state.geometry,
                    packet_source,
                    legacy_mode_enabled,
                )
        if "plasma" in kwargs:
            plasma = kwargs["plasma"]
        else:
            plasma = assemble_plasma(
                config,
                simulation_state,
                atom_data=atom_data,
            )
        if "transport" in kwargs:
            if packet_source is not None:
                raise ConfigurationError(
                    "Cannot specify packet_source and transport at the same time."
                )
            transport = kwargs["transport"]
        else:
            transport = MonteCarloTransportSolver.from_config(
                config,
                packet_source=simulation_state.packet_source,
                enable_virtual_packet_logging=virtual_packet_logging,
            )

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

        return cls(
            iterations=config.montecarlo.iterations,
            simulation_state=simulation_state,
            plasma=plasma,
            transport=transport,
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
        )
