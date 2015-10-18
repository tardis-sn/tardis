import logging
import time
from pandas import HDFStore
import os

import numpy as np

from tardis.model import Radial1DModel
from tardis.montecarlo.base import MontecarloRunner

# Adding logging support
logger = logging.getLogger(__name__)


class Simulation(object):
    def __init__(self, tardis_config):
        self.tardis_config = tardis_config
        self.model = Radial1DModel(tardis_config)
        self.runner = MontecarloRunner(self.tardis_config.montecarlo.seed)

    def do_montecarlo(self, no_of_packets, no_of_virtual_packets):
        if update_radiation_field:
            t_inner_new = self.update_radiationfield()
        else:
            t_inner_new = self.t_inner

        self.calculate_j_blues(init_detailed_j_blues=initialize_j_blues)
        self.update_plasmas(initialize_nlte=initialize_nlte)

        self.t_inner = t_inner_new


        self.montecarlo_virtual_luminosity = np.zeros_like(self.spectrum.frequency.value)

        self.runner.run(self, no_of_virtual_packets=no_of_virtual_packets,
                        nthreads=self.tardis_config.montecarlo.nthreads) #self = model


        (montecarlo_nu, montecarlo_energies, self.j_estimators,
         self.nubar_estimators, last_line_interaction_in_id,
         last_line_interaction_out_id, self.last_interaction_type,
         self.last_line_interaction_shell_id) = self.runner.legacy_return()

        if np.sum(montecarlo_energies < 0) == len(montecarlo_energies):
            logger.critical("No r-packet escaped through the outer boundary.")

        self.montecarlo_nu = self.runner.output_nu
        self.montecarlo_luminosity = self.runner.packet_luminosity



        montecarlo_reabsorbed_luminosity = np.histogram(
            self.runner.reabsorbed_packet_nu,
            weights=self.runner.reabsorbed_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s



        montecarlo_emitted_luminosity = np.histogram(
            self.runner.emitted_packet_nu,
            weights=self.runner.emitted_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s



        self.spectrum.update_luminosity(montecarlo_emitted_luminosity)
        self.spectrum_reabsorbed.update_luminosity(montecarlo_reabsorbed_luminosity)


        if no_of_virtual_packets > 0:
            self.montecarlo_virtual_luminosity = self.montecarlo_virtual_luminosity \
                                                 * 1 * u.erg / self.time_of_simulation
            self.spectrum_virtual.update_luminosity(self.montecarlo_virtual_luminosity)



        self.last_line_interaction_in_id = self.atom_data.lines_index.index.values[last_line_interaction_in_id]
        self.last_line_interaction_in_id = self.last_line_interaction_in_id[last_line_interaction_in_id != -1]
        self.last_line_interaction_out_id = self.atom_data.lines_index.index.values[last_line_interaction_out_id]
        self.last_line_interaction_out_id = self.last_line_interaction_out_id[last_line_interaction_out_id != -1]
        self.last_line_interaction_angstrom = self.montecarlo_nu[last_line_interaction_in_id != -1].to('angstrom',
                                                                                                       u.spectral())


        self.iterations_executed += 1
        self.iterations_remaining -= 1

        if self.gui is not None:
            self.gui.update_data(self)
            self.gui.show()


    def simulate(self, update_radiation_field=True, enable_virtual=False, initialize_j_blues=False,
                 initialize_nlte=False):
        """
        Run a simulation
        """





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










