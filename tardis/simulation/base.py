import logging
import time
import itertools

import pandas as pd

import numpy as np

from astropy import units as u

from tardis.montecarlo.base import MontecarloRunner

# Adding logging support
logger = logging.getLogger(__name__)


class Simulation(object):

    converged = False

    def __init__(self, tardis_config):
        self.tardis_config = tardis_config
        self.runner = MontecarloRunner(self.tardis_config.montecarlo.seed,
                                       tardis_config.spectrum.frequency)
        t_inner_lock_cycle = [False] * (tardis_config.montecarlo.
                                        convergence_strategy.
                                        lock_t_inner_cycles)
        t_inner_lock_cycle[0] = True
        self.t_inner_update = itertools.cycle(t_inner_lock_cycle)


    def run_single_montecarlo(self, model, no_of_packets,
                              no_of_virtual_packets=0):
        """
        Will do a single TARDIS iteration with the given model
        Parameters
        ----------
        model: ~tardis.model.Radial1DModel
        no_of_packet: ~int
        no_of_virtual_packets: ~int
            default is 0 and switches of the virtual packet mode. Recommended
            is 3.

        Returns
        -------
            : None

        """
        self.runner.run(model, no_of_packets,
                        no_of_virtual_packets=no_of_virtual_packets,
                        nthreads=self.tardis_config.montecarlo.nthreads)


        (montecarlo_nu, montecarlo_energies, self.j_estimators,
         self.nubar_estimators, last_line_interaction_in_id,
         last_line_interaction_out_id, self.last_interaction_type,
         self.last_line_interaction_shell_id) = self.runner.legacy_return()

        if np.sum(montecarlo_energies < 0) == len(montecarlo_energies):
            logger.critical("No r-packet escaped through the outer boundary.")



    def calculate_emitted_luminosity(self):
        """

        Returns
        -------

        """
        return self.runner.calculate_emitted_luminosity(
            self.tardis_config.supernova.luminosity_nu_start,
            self.tardis_config.supernova.luminosity_nu_end)

    def calculate_reabsorbed_luminosity(self):
        return self.runner.calculate_reabsorbed_luminosity(
            self.tardis_config.supernova.luminosity_nu_start,
            self.tardis_config.supernova.luminosity_nu_end)


    def estimate_t_inner(self, input_t_inner, luminosity_requested,
                         t_inner_update_exponent=-0.5):
        emitted_luminosity = self.calculate_emitted_luminosity()

        luminosity_ratios = (
            (emitted_luminosity / luminosity_requested).to(1).value)

        return input_t_inner * luminosity_ratios ** t_inner_update_exponent

    def get_convergence_status(self, t_rad, w, t_inner, estimated_t_rad, estimated_w,
                               estimated_t_inner):
        convergence_section = self.tardis_config.montecarlo.convergence_strategy
        no_of_shells = self.tardis_config.structure.no_of_shells

        convergence_t_rad = (abs(t_rad - estimated_t_rad) /
                             estimated_t_rad).value
        convergence_w = (abs(w - estimated_w) / estimated_w)
        convergence_t_inner = (abs(t_inner - estimated_t_inner) /
                               estimated_t_inner).value

        if convergence_section.type == 'specific':
            fraction_t_rad_converged = (
                np.count_nonzero(
                    convergence_t_rad < convergence_section.t_rad.threshold)
                / no_of_shells)

            t_rad_converged = (
                fraction_t_rad_converged > convergence_section.t_rad.threshold)

            fraction_w_converged = (
                np.count_nonzero(
                    convergence_w < convergence_section.w.threshold)
                / no_of_shells)

            w_converged = (
                fraction_w_converged > convergence_section.w.threshold)

            t_inner_converged = (
                convergence_t_inner < convergence_section.t_inner.threshold)

            if np.all([t_rad_converged, w_converged, t_inner_converged]):
                return True
            else:
                return False

        else:
            return False


    def log_run_results(self, emitted_luminosity, absorbed_luminosity):
            logger.info("Luminosity emitted = {0:.5e} "
                    "Luminosity absorbed = {1:.5e} "
                    "Luminosity requested = {2:.5e}".format(
                emitted_luminosity, absorbed_luminosity,
                self.tardis_config.supernova.luminosity_requested))


    def log_plasma_state(self, t_rad, w, t_inner, next_t_rad, next_w,
                         next_t_inner, log_sampling=5):
        """
        Logging the change of the plasma state

        Parameters
        ----------
        t_rad: ~astropy.units.Quanity
            current t_rad
        w: ~astropy.units.Quanity
            current w
        next_t_rad: ~astropy.units.Quanity
            next t_rad
        next_w: ~astropy.units.Quanity
            next_w
        log_sampling: ~int
            the n-th shells to be plotted

        Returns
        -------

        """

        plasma_state_log = pd.DataFrame(index=np.arange(len(t_rad)),
                                           columns=['t_rad', 'next_t_rad',
                                                    'w', 'next_w'])
        plasma_state_log['t_rad'] = t_rad
        plasma_state_log['next_t_rad'] = next_t_rad
        plasma_state_log['w'] = w
        plasma_state_log['next_w'] = next_w

        plasma_state_log.index.name = 'Shell'

        plasma_state_log = str(plasma_state_log[::log_sampling])

        plasma_state_log = ''.join(['\t%s\n' % item for item in
                                    plasma_state_log.split('\n')])

        logger.info('Plasma stratification:\n%s\n', plasma_state_log)
        logger.info('t_inner {0:.3f} -- next t_inner {1:.3f}'.format(
            t_inner, next_t_inner))


    @staticmethod
    def damped_converge(value, estimated_value, damping_factor):
        return value + damping_factor * (estimated_value - value)


    def calculate_next_plasma_state(self, t_rad, w, t_inner,
                                    estimated_w, estimated_t_rad,
                                    estimated_t_inner):

        convergence_strategy = (
            self.tardis_config.montecarlo.convergence_strategy)

        if (convergence_strategy.type == 'damped'
            or convergence_strategy.type == 'specific'):

            next_t_rad = self.damped_converge(
                t_rad, estimated_t_rad,
                convergence_strategy.t_rad.damping_constant)
            next_w = self.damped_converge(
                w, estimated_w, convergence_strategy.w.damping_constant)
            next_t_inner = self.damped_converge(
                t_inner, estimated_t_inner,
                convergence_strategy.t_inner.damping_constant)

            return next_t_rad, next_w, next_t_inner

        else:
            raise ValueError('Convergence strategy type is '
                             'neither damped nor specific '
                             '- input is {0}'.format(convergence_strategy.type))

    def legacy_run_simulation(self, model):
        start_time = time.time()

        iterations_remaining = self.tardis_config.montecarlo.iterations
        iterations_max_requested = self.tardis_config.montecarlo.iterations
        iterations_executed = 0
        converged = False

        while iterations_remaining > 1:
            logger.info('Remaining run %d', iterations_remaining)
            self.run_single_montecarlo(
                model, self.tardis_config.montecarlo.no_of_packets)
            self.log_run_results(self.calculate_emitted_luminosity(),
                                 self.calculate_reabsorbed_luminosity())
            iterations_executed += 1
            iterations_remaining -= 1

            estimated_t_rad, estimated_w = (
                self.runner.calculate_radiationfield_properties())
            estimated_t_inner = self.estimate_t_inner(
                model.t_inner,
                self.tardis_config.supernova.luminosity_requested)

            converged = self.get_convergence_status(
                model.t_rads, model.ws, model.t_inner, estimated_t_rad,
                estimated_w, estimated_t_inner)

            next_t_rad, next_w, next_t_inner = self.calculate_next_plasma_state(
                model.t_rads, model.ws, model.t_inner,
                estimated_w, estimated_t_rad, estimated_t_inner)

            self.log_plasma_state(model.t_rads, model.ws, model.t_inner,
                                  next_t_rad, next_w, next_t_inner)
            model.t_rads = next_t_rad
            model.ws = next_w
            model.t_inner = next_t_inner
            model.j_blue_estimators = self.runner.j_blue_estimator

            model.calculate_j_blues(init_detailed_j_blues=False)
            model.update_plasmas(initialize_nlte=False)


            # if switching into the hold iterations mode or out back to the normal one
            # if it is in either of these modes already it will just stay there
            if converged and not self.converged:
                self.converged = True
                iterations_remaining = (
                    convergence_section.global_convergence_parameters.
                        hold_iterations_wrong)
            elif not converged and self.converged:
                # UMN Warning: the following two iterations attributes of the Simulation object don't exist
                self.iterations_remaining = self.iterations_max_requested - self.iterations_executed
                self.converged = False
            else:
                # either it is converged and the status of the simulation is
                # converged OR it is not converged and the status of the
                # simulation is not converged - Do nothing.
                pass

            if converged:
                convergence_section = (
                    self.tardis_config.montecarlo.convergence_strategy)
                iterations_remaining = (
                    convergence_section.global_convergence_parameters.
                        hold_iterations)

        #Finished second to last loop running one more time
        logger.info('Doing last run')
        if self.tardis_config.montecarlo.last_no_of_packets is not None:
            no_of_packets = self.tardis_config.montecarlo.last_no_of_packets
        else:
            no_of_packets = self.tardis_config.montecarlo.no_of_packets

        no_of_virtual_packets = (
            self.tardis_config.montecarlo.no_of_virtual_packets)

        self.run_single_montecarlo(model, no_of_packets, no_of_virtual_packets)

        self.legacy_update_spectrum(model, no_of_virtual_packets)
        self.legacy_set_final_model_properties(model)

        #the following instructions, passing down information to the model are
        #required for the gui
        model.no_of_packets = no_of_packets
        model.no_of_virtual_packets = no_of_virtual_packets
        model.converged = converged
        model.iterations_executed = iterations_executed
        model.iterations_max_requested = iterations_max_requested

        logger.info("Finished in {0:d} iterations and took {1:.2f} s".format(
            iterations_executed, time.time()-start_time))


    def legacy_update_spectrum(self, model, no_of_virtual_packets):
        montecarlo_reabsorbed_luminosity = np.histogram(
            self.runner.reabsorbed_packet_nu,
            weights=self.runner.reabsorbed_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s

        montecarlo_emitted_luminosity = np.histogram(
            self.runner.emitted_packet_nu,
            weights=self.runner.emitted_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s

        model.spectrum.update_luminosity(montecarlo_emitted_luminosity)
        model.spectrum_reabsorbed.update_luminosity(montecarlo_reabsorbed_luminosity)

        if no_of_virtual_packets > 0:
            model.montecarlo_virtual_luminosity = (
                self.runner.legacy_montecarlo_virtual_luminosity *
                1 * u.erg / model.time_of_simulation)[:-1]
            model.spectrum_virtual.update_luminosity(
                model.montecarlo_virtual_luminosity)

    def legacy_set_final_model_properties(self, model):
        """Sets additional model properties to be compatible with old model design

        The runner object is given to the model and other packet diagnostics are set.

        Parameters
        ----------
        model: ~tardis.model.Radial1DModel

        Returns
        -------
            : None

        """

        #pass the runner to the model
        model.runner = self.runner
        #TODO: pass packet diagnostic arrays
        (montecarlo_nu, montecarlo_energies, model.j_estimators,
                model.nubar_estimators, last_line_interaction_in_id,
                last_line_interaction_out_id, model.last_interaction_type,
                model.last_line_interaction_shell_id) = model.runner.legacy_return()

        model.montecarlo_nu = self.runner.output_nu
        model.montecarlo_luminosity = self.runner.packet_luminosity


        model.last_line_interaction_in_id = model.atom_data.lines_index.index.values[last_line_interaction_in_id]
        model.last_line_interaction_in_id = model.last_line_interaction_in_id[last_line_interaction_in_id != -1]
        model.last_line_interaction_out_id = model.atom_data.lines_index.index.values[last_line_interaction_out_id]
        model.last_line_interaction_out_id = model.last_line_interaction_out_id[last_line_interaction_out_id != -1]
        model.last_line_interaction_angstrom = model.montecarlo_nu[last_line_interaction_in_id != -1].to('angstrom',
                                                                                                       u.spectral())
        # required for gui
        model.current_no_of_packets = model.tardis_config.montecarlo.no_of_packets

def run_radial1d(radial1d_model, history_fname=None):
    if history_fname is not None:
        raise ValueError('This functionality is currently not supported')

    simulation = Simulation(radial1d_model.tardis_config)
    simulation.legacy_run_simulation(radial1d_model)










