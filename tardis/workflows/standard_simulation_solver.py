import logging

import numpy as np
from astropy import units as u

from tardis.io.atom_data.base import AtomData
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.formal_integral import FormalIntegrator

# logging support
logger = logging.getLogger(__name__)


class StandardSimulationSolver:
    def __init__(self, convergence_strategy, atom_data_path):
        self.atom_data_path = atom_data_path
        self.simulation_state = None
        self.atom_data = None
        self.spectrum_solver = None
        self.transport_solver = None
        self.plasma_solver = None
        self.luminosity_requested = 0 * u.erg / u.s
        self.integrated_spectrum_settings = None

        # Convergence
        self.convergence_strategy = convergence_strategy
        self.consecutive_converges_count = 0
        self.converged = False
        self.total_iterations = 1
        self.completed_iterations = 0

        # Convergence solvers
        self.t_radiative_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_radiative
        )
        self.dilution_factor_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.dilution_factor
        )
        self.t_inner_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_inner
        )

    def _get_convergence_status(
        self,
        t_radiative,
        dilution_factor,
        t_inner,
        estimated_t_radiative,
        estimated_dilution_factor,
        estimated_t_inner,
    ):
        t_radiative_converged = (
            self.t_radiative_convergence_solver.get_convergence_status(
                t_radiative.value,
                estimated_t_radiative.value,
                self.simulation_state.no_of_shells,
            )
        )

        dilution_factor_converged = (
            self.dilution_factor_convergence_solver.get_convergence_status(
                dilution_factor,
                estimated_dilution_factor,
                self.simulation_state.no_of_shells,
            )
        )

        t_inner_converged = (
            self.t_inner_convergence_solver.get_convergence_status(
                t_inner.value,
                estimated_t_inner.value,
                1,
            )
        )

        if np.all(
            [
                t_radiative_converged,
                dilution_factor_converged,
                t_inner_converged,
            ]
        ):
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

        self.consecutive_converges_count = 0
        return False

    def _get_atom_data(self):
        try:
            self.atom_data = AtomData.from_hdf(self.atom_data_path)
        except TypeError:
            logger.debug("Atom Data Cannot be Read from HDF.")

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

    def get_convergence_estimates(self, emitted_luminosity):
        (
            estimated_t_radiative,
            estimated_dilution_factor,
        ) = self.transport_solver.transport_state.calculate_radiationfield_properties()

        estimated_t_inner = self.estimate_t_inner(
            self.simulation_state.t_inner,
            self.luminosity_requested,
            emitted_luminosity,
            t_inner_update_exponent=self.convergence_strategy.t_inner_update_exponent,
        )
        return (
            estimated_t_radiative,
            estimated_dilution_factor,
            estimated_t_inner,
        )

    def check_convergence(
        self,
        estimated_t_radiative,
        estimated_dilution_factor,
        estimated_t_inner,
    ):
        converged = self._get_convergence_status(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_radiative,
            estimated_dilution_factor,
            estimated_t_inner,
        )

        return converged

    def solve_plasma(
        self,
        estimated_t_radiative,
        estimated_dilution_factor,
        estimated_t_inner,
    ):
        next_t_radiative = self.t_rad_convergence_solver.converge(
            self.simulation_state.t_radiative,
            estimated_t_radiative,
        )
        next_dilution_factor = self.dilution_factor_convergence_solver.converge(
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

        self.simulation_state.t_radiative = next_t_radiative
        self.simulation_state.dilution_factor = next_dilution_factor
        self.simulation_state.blackbody_packet_source.temperature = next_t_inner

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

        self.plasma_solver.update(**update_properties)

    def solve_montecarlo(self, no_of_real_packets, no_of_virtual_packets):
        transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport.run(
            transport_state,
            time_explosion=self.simulation_state.time_explosion,
            iteration=self.completed_iterations,
            total_iterations=self.total_iterations,
            show_progress_bars=self.show_progress_bars,
        )

        return transport_state, virtual_packet_energies

    def solve_spectrum(
        self,
        transport_state,
        virtual_packet_energies=None,
        integrated_spectrum_settings=None,
    ):
        # Set up spectrum solver
        self.spectrum_solver.transport_state = transport_state
        if virtual_packet_energies is not None:
            self.spectrum_solver._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_energies
            )

        if integrated_spectrum_settings is not None:
            # Set up spectrum solver integrator
            self.spectrum_solver.integrator_settings = (
                integrated_spectrum_settings
            )
            self.spectrum_solver._integrator = FormalIntegrator(
                self.simulation_state, self.plasma, self.transport
            )

    def calculate_emitted_luminosity(self, transport_state):
        self.spectrum_solver.transport_state = transport_state

        output_energy = (
            self.transport.transport_state.packet_collection.output_energies
        )
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        emitted_luminosity = self.spectrum_solver.calculate_emitted_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end
        )
        return emitted_luminosity

    def solve(self):
        converged = False
        while self.completed_iterations < self.total_iterations - 1:
            transport_state, virtual_packet_energies = self.solve_montecarlo()

            emitted_luminosity = self.calculate_emitted_luminosity(
                transport_state
            )

            (
                estimated_t_radiative,
                estimated_dilution_factor,
                estimated_t_inner,
            ) = self.get_convergence_estimates(emitted_luminosity)

            self.solve_plasma(
                estimated_t_radiative,
                estimated_dilution_factor,
                estimated_t_inner,
            )

            converged = self.check_convergence(
                estimated_t_radiative,
                estimated_dilution_factor,
                estimated_t_inner,
            )
            self.completed_iterations += 1

            if converged and self.convergence_strategy.stop_if_converged:
                break

        transport_state, virtual_packet_energies = self.solve_montecarlo()
        self.solve_spectrum(
            transport_state,
            virtual_packet_energies,
            self.integrated_spectrum_settings,
        )
