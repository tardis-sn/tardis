import logging
import time
import os
import itertools

from pandas import HDFStore

import numpy as np

from tardis.model import Radial1DModel
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

        Parameters
        ----------
        no_of_packets
        no_of_virtual_packets

        Returns
        -------

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

    def estimate_new_t_inner(self, input_t_inner, luminosity_requested,
                             t_inner_update_exponent=0.5):
        emitted_luminosity = self.runner.calculate_emitted_luminosity(
            self.tardis_config.supernova.luminosity_nu_start,
            self.tardis_config.supernova.luminosity_nu_end)

        luminosity_ratios =(
            (emitted_luminosity / luminosity_requested).to(1).value)

        return input_t_inner * luminosity_ratios ** t_inner_update_exponent

    def get_convergence_status(self, t_rad, w, t_inner, new_t_rad, new_w,
                               new_t_inner):
        convergence_section = self.tardis_config.montecarlo.convergence_strategy
        no_of_shells = self.tardis_config.structure.no_of_shells

        convergence_t_rad = (abs(t_rad - new_t_rad) / new_t_rad).value
        convergence_w = (abs(w - new_w) / new_w)
        convergence_t_inner = (abs(t_inner - new_t_inner) / new_t_inner).value

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





    def legacy_run_simulation(self, model):
        start_time = time.time()

        iterations_remaining = self.tardis_config.montecarlo.iterations
        iterations_executed = 0

        while iterations_remaining > 1:
            logger.info('Remaining run %d', iterations_remaining)
            self.run_single_montecarlo(
                model, self.tardis_config.montecarlo.no_of_packets)
            iterations_executed += 1

            new_t_rad, new_w = self.runner.calculate_radiationfield_properties()
            new_t_inner = self.estimate_new_t_inner()

            converged = self.get_convergence_status(
                model.t_rads, model.ws, model.t_inner, new_t_rad, new_w,
                new_t_inner)

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

        logger.info("Finished in {0:d} iterations and took {1:.2f} s".format(
            iterations_executed, time.time()-start_time))

    def update_radiationfield(self, log_sampling=5):
        """
        Updating radiation field
        """
        convergence_section = self.tardis_config.montecarlo.convergence_strategy
        updated_t_rads, updated_ws = (
            self.runner.calculate_radiationfield_properties())
        old_t_rads = self.t_rads.copy()
        old_ws = self.ws.copy()
        old_t_inner = self.t_inner
        luminosity_wavelength_filter = (self.montecarlo_nu > self.tardis_config.supernova.luminosity_nu_start) & \
                            (self.montecarlo_nu < self.tardis_config.supernova.luminosity_nu_end)
        emitted_filter = self.montecarlo_luminosity.value >= 0
        emitted_luminosity = np.sum(self.montecarlo_luminosity.value[emitted_filter & luminosity_wavelength_filter]) \
                             * self.montecarlo_luminosity.unit

        absorbed_luminosity = -np.sum(self.montecarlo_luminosity.value[~emitted_filter & luminosity_wavelength_filter]) \
                              * self.montecarlo_luminosity.unit
        updated_t_inner = self.t_inner \
                          * (emitted_luminosity / self.tardis_config.supernova.luminosity_requested).to(1).value \
                            ** convergence_section.t_inner_update_exponent

        if convergence_section.type == 'damped' or convergence_section.type == 'specific':
            self.t_rads += convergence_section.t_rad.damping_constant * (updated_t_rads - self.t_rads)
            self.ws += convergence_section.w.damping_constant * (updated_ws - self.ws)
            if self.t_inner_update.next():
                t_inner_new = self.t_inner + convergence_section.t_inner.damping_constant * (updated_t_inner - self.t_inner)
            else:
                t_inner_new = self.t_inner


        if convergence_section.type == 'specific':

            if t_rad_converged and t_inner_converged and w_converged:
                if not self.converged:
                    self.converged = True


            else:
                if self.converged:
                    self.iterations_remaining = self.iterations_max_requested - self.iterations_executed
                    self.converged = False

        self.temperature_logging = pd.DataFrame(
            {'t_rads': old_t_rads.value, 'updated_t_rads': updated_t_rads.value,
             'converged_t_rads': convergence_t_rads, 'new_trads': self.t_rads.value, 'ws': old_ws,
             'updated_ws': updated_ws, 'converged_ws': convergence_ws,
             'new_ws': self.ws})

        self.temperature_logging.index.name = 'Shell'

        temperature_logging = str(self.temperature_logging[::log_sampling])

        temperature_logging = ''.join(['\t%s\n' % item for item in temperature_logging.split('\n')])

        logger.info('Plasma stratification:\n%s\n', temperature_logging)
        logger.info("Luminosity emitted = %.5e Luminosity absorbed = %.5e Luminosity requested = %.5e",
                    emitted_luminosity.value, absorbed_luminosity.value,
                    self.tardis_config.supernova.luminosity_requested.value)
        logger.info('Calculating new t_inner = %.3f', updated_t_inner.value)

        return t_inner_new


def run_radial1d(radial1d_model, history_fname=None):
    if history_fname:
        if os.path.exists(history_fname):
            logger.warn('History file %s exists - it will be overwritten', history_fname)
            os.system('rm %s' % history_fname)
        history_buffer = HDFStore(history_fname)
        radial1d_model.atom_data.lines.to_hdf(history_buffer, 'atom_data/lines')
        radial1d_model.atom_data.levels.to_hdf(history_buffer, 'atom_data/levels')


    start_time = time.time()
    initialize_j_blues = True
    initialize_nlte = True
    update_radiation_field = False
    while radial1d_model.iterations_remaining > 1:
        logger.info('Remaining run %d', radial1d_model.iterations_remaining)
        radial1d_model.simulate(update_radiation_field=update_radiation_field, enable_virtual=False, initialize_nlte=initialize_nlte,
                                initialize_j_blues=initialize_j_blues)
        initialize_j_blues=False
        initialize_nlte=False
        update_radiation_field = True

        if history_fname:
            radial1d_model.to_hdf5(history_buffer, path='model%03d' % radial1d_model.iterations_executed, close_h5=False)

    #Finished second to last loop running one more time
    logger.info('Doing last run')
    if radial1d_model.tardis_config.montecarlo.last_no_of_packets is not None:
        radial1d_model.current_no_of_packets = radial1d_model.tardis_config.montecarlo.last_no_of_packets

    radial1d_model.simulate(enable_virtual=True, update_radiation_field=update_radiation_field, initialize_nlte=initialize_nlte,
                            initialize_j_blues=initialize_j_blues)

    if history_fname:
        radial1d_model.to_hdf5(history_buffer, path='model%03d' % radial1d_model.iterations_executed)



    logger.info("Finished in %d iterations and took %.2f s", radial1d_model.iterations_executed, time.time()-start_time)










