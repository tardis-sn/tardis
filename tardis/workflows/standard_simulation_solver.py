import logging
from pathlib import Path

import numpy as np
from astropy import units as u

from tardis import constants as const
from tardis.io.atom_data.base import AtomData
from tardis.model import SimulationState
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral import FormalIntegrator
from tardis.transport.montecarlo.base import MonteCarloTransportSolver

# logging support
logger = logging.getLogger(__name__)


class StandardSimulationSolver:
    def __init__(self, configuration):
        # Convergence
        self.consecutive_converges_count = 0
        self.converged = False
        self.completed_iterations = 0
        self.luminosity_requested = (
            configuration.supernova.luminosity_requested.cgs
        )

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
            enable_virtual_packet_logging=False,
        )

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

    def _get_atom_data(self, configuration):
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
        (
            estimated_t_radiative,
            estimated_dilution_factor,
        ) = self.transport_solver.transport_state.calculate_radiationfield_properties()

        self.initialize_spectrum_solver(
            transport_state,
            None,
        )

        emitted_luminosity = self.spectrum_solver.calculate_emitted_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end
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
        }

    def check_convergence(
        self,
        estimated_values,
    ):
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

    def solve_simulation_state(
        self,
        estimated_values,
    ):
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

    def solve_plasma(
        self,
        transport_state,
    ):
        update_properties = dict(
            t_rad=self.simulation_state.t_radiative,
            w=self.simulation_state.dilution_factor,
        )
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.
        if "j_blue_estimator" in self.plasma_solver.outputs_dict:
            update_properties.update(
                t_inner=self.simulation_state.blackbody_packet_source.temperature,
                j_blue_estimator=transport_state.radfield_mc_estimators.j_blue_estimator,
            )

        self.plasma_solver.update(**update_properties)

    def solve_montecarlo(self, no_of_real_packets, no_of_virtual_packets=0):
        transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport_solver.run(
            transport_state,
            time_explosion=self.simulation_state.time_explosion,
            iteration=self.completed_iterations,
            total_iterations=self.total_iterations,
            show_progress_bars=False,
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
        # Set up spectrum solver
        self.spectrum_solver.transport_state = transport_state

        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_energies
            )

        if self.integrated_spectrum_settings is not None:
            # Set up spectrum solver integrator
            self.spectrum_solver.integrator_settings = (
                self.integrated_spectrum_settings
            )
            self.spectrum_solver._integrator = FormalIntegrator(
                self.simulation_state, self.plasma_solver, self.transport_solver
            )

    def solve(self):
        converged = False
        while self.completed_iterations < self.total_iterations - 1:
            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.real_packet_count
            )

            estimated_values = self.get_convergence_estimates(transport_state)

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(transport_state)

            converged = self.check_convergence(estimated_values)
            self.completed_iterations += 1

            if converged and self.convergence_strategy.stop_if_converged:
                break

        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count, self.virtual_packet_count
        )
        self.initialize_spectrum_solver(
            transport_state,
            virtual_packet_energies,
        )
