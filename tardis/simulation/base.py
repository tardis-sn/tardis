import time
import logging
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.montecarlo.base import MontecarloRunner
from tardis.model.base import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma

# Adding logging support
logger = logging.getLogger(__name__)


class Simulation(object):
    def __init__(self, iterations, hold_iterations, model, plasma, runner,
                 no_of_packets, no_of_virtual_packets, luminosity_nu_start,
                 luminosity_nu_end, last_no_of_packets,
                 luminosity_requested, convergence_strategy,
                 nthreads):
        self.converged = False
        self.iterations = iterations
        self.hold_iterations = hold_iterations
        self.iterations_executed = 0
        self.model = model
        self.plasma = plasma
        self.runner = runner
        self.no_of_packets = no_of_packets
        self.last_no_of_packets = last_no_of_packets
        self.no_of_virtual_packets = no_of_virtual_packets
        self.luminosity_nu_start = luminosity_nu_start
        self.luminosity_nu_end = luminosity_nu_end
        self.luminosity_requested = luminosity_requested
        self.nthreads = nthreads
        if convergence_strategy.type in ('damped', 'specific'):
            self.convergence_strategy = convergence_strategy
        else:
            raise ValueError('Convergence strategy type is '
                             'neither damped nor specific '
                             '- input is {0}'.format(convergence_strategy.type))

    def estimate_t_inner(self, input_t_inner, luminosity_requested,
                         t_inner_update_exponent=-0.5):
        emitted_luminosity = self.runner.calculate_emitted_luminosity(
                                self.luminosity_nu_start,
                                self.luminosity_nu_end)

        luminosity_ratios = (
            (emitted_luminosity / luminosity_requested).to(1).value)

        return input_t_inner * luminosity_ratios ** t_inner_update_exponent

    @staticmethod
    def damped_converge(value, estimated_value, damping_factor):
        # FIXME: Should convergence strategy have its own class containing this
        # as a method
        return value + damping_factor * (estimated_value - value)

    def _get_convergence_status(self, t_rad, w, t_inner, estimated_t_rad,
                                estimated_w, estimated_t_inner):
        # FIXME: Move the convergence checking in its own class.
        no_of_shells = self.model.no_of_shells

        convergence_t_rad = (abs(t_rad - estimated_t_rad) /
                             estimated_t_rad).value
        convergence_w = (abs(w - estimated_w) / estimated_w)
        convergence_t_inner = (abs(t_inner - estimated_t_inner) /
                               estimated_t_inner).value

        if self.convergence_strategy.type == 'specific':
            fraction_t_rad_converged = (
                np.count_nonzero(
                    convergence_t_rad < self.convergence_strategy.t_rad.threshold)
                / no_of_shells)

            t_rad_converged = (
                fraction_t_rad_converged > self.convergence_strategy.t_rad.threshold)

            fraction_w_converged = (
                np.count_nonzero(
                    convergence_w < self.convergence_strategy.w.threshold)
                / no_of_shells)

            w_converged = (
                fraction_w_converged > self.convergence_strategy.w.threshold)

            t_inner_converged = (
                convergence_t_inner < self.convergence_strategy.t_inner.threshold)

            if np.all([t_rad_converged, w_converged, t_inner_converged]):
                return True
            else:
                return False

        else:
            return False

    def advance_state(self):
        """
        Advances the state of the model and the plasma for the next
        iteration of the simulation. Returns True if the convergence criteria
        are met, else False.

        Returns
        -------
            converged : ~bool
        """
        estimated_t_rad, estimated_w = (
            self.runner.calculate_radiationfield_properties())
        estimated_t_inner = self.estimate_t_inner(
            self.model.t_inner, self.luminosity_requested)

        converged = self._get_convergence_status(self.model.t_rad,
                                                 self.model.w,
                                                 self.model.t_inner,
                                                 estimated_t_rad,
                                                 estimated_w,
                                                 estimated_t_inner)

        # calculate_next_plasma_state equivalent
        # FIXME: Should convergence strategy have its own class?
        next_t_rad = self.damped_converge(
            self.model.t_rad, estimated_t_rad,
            self.convergence_strategy.t_rad.damping_constant)
        next_w = self.damped_converge(
            self.model.w, estimated_w, self.convergence_strategy.w.damping_constant)
        next_t_inner = self.damped_converge(
            self.model.t_inner, estimated_t_inner,
            self.convergence_strategy.t_inner.damping_constant)

        self.log_plasma_state(self.model.t_rad, self.model.w,
                              self.model.t_inner, next_t_rad, next_w,
                              next_t_inner)
        self.model.t_rad = next_t_rad
        self.model.w = next_w
        self.model.t_inner = next_t_inner

        # model.calculate_j_blues() equivalent
        # model.update_plasmas() equivalent
        # TODO: if nlte_config is not None and nlte_config.species:
        #          self.store_previous_properties()
        self.plasma.update(t_rad=self.model.t_rad, w=self.model.w)

        return converged

    def iterate(self, no_of_packets, no_of_virtual_packets=0, last_run=False):
        logger.info('Starting iteration {0:d}/{1:d}'.format(
                    self.iterations_executed + 1, self.iterations))
        self.runner.run(self.model, self.plasma, no_of_packets,
                        no_of_virtual_packets=no_of_virtual_packets,
                        nthreads=self.nthreads, last_run=last_run)
        output_energy = self.runner.output_energy
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        emitted_luminosity = self.runner.calculate_emitted_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end)
        reabsorbed_luminosity = self.runner.calculate_reabsorbed_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end)
        self.log_run_results(emitted_luminosity,
                             reabsorbed_luminosity)
        self.iterations_executed += 1

    def run(self):
        start_time = time.time()
        times_converged = 0
        # If an iteration has converged, require hold_iterations more iterations
        # to converge before we conclude that the Simulation is converged.
        while (self.iterations_executed < self.iterations - 1 and
                        times_converged <= self.hold_iterations):
            self.iterate(self.no_of_packets)
            converged = self.advance_state()
            if converged:
                times_converged += 1
                logger.info("Iteration converged {0:d}/{1:d} consecutive "
                            "times.".format(times_converged,
                                            self.hold_iterations + 1))
            else:
                times_converged = 0
        self.converged = times_converged == self.hold_iterations + 1
        # Last iteration
        self.iterate(self.last_no_of_packets, self.no_of_virtual_packets, True)
        self.runner.legacy_update_spectrum(self.no_of_virtual_packets)

        logger.info("Simulation finished in {0:d} iterations "
                    "and took {1:.2f} s".format(
                        self.iterations_executed, time.time() - start_time))

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

    def log_run_results(self, emitted_luminosity, absorbed_luminosity):
        logger.info("Luminosity emitted = {0:.5e} "
                    "Luminosity absorbed = {1:.5e} "
                    "Luminosity requested = {2:.5e}".format(
            emitted_luminosity, absorbed_luminosity,
            self.luminosity_requested))

    def to_hdf(self, path_or_buf, path='', plasma_properties=None):
        """
        Store the simulation to an HDF structure.

        Parameters
        ----------
        path_or_buf
            Path or buffer to the HDF store
        path : str
            Path inside the HDF store to store the simulation
        plasma_properties
            `None` or a `PlasmaPropertyCollection` which will
            be passed as the collection argument to the
            plasma.to_hdf method.
        Returns
        -------
        None
        """
        self.runner.to_hdf(path_or_buf, path)
        self.model.to_hdf(path_or_buf, path)
        self.plasma.to_hdf(path_or_buf, path, plasma_properties)

    @classmethod
    def from_config(cls, config):
        model = Radial1DModel.from_config(config)
        plasma = assemble_plasma(config, model)
        runner = MontecarloRunner.from_config(config)

        luminosity_nu_start = config.supernova.luminosity_wavelength_end.to(
                u.Hz, u.spectral())

        try:
            luminosity_nu_end = config.supernova.luminosity_wavelength_start.to(
                u.Hz, u.spectral())
        except ZeroDivisionError:
            luminosity_nu_end = np.inf * u.Hz

        last_no_of_packets = config.montecarlo.last_no_of_packets
        if last_no_of_packets is None or last_no_of_packets < 0:
            last_no_of_packets =  config.montecarlo.no_of_packets
        last_no_of_packets = int(last_no_of_packets)

        return cls(iterations=config.montecarlo.iterations,
                   hold_iterations=config.montecarlo.convergence_strategy.hold_iterations,
                   model=model,
                   plasma=plasma,
                   runner=runner,
                   no_of_packets=int(config.montecarlo.no_of_packets),
                   no_of_virtual_packets=int(
                       config.montecarlo.no_of_virtual_packets),
                   luminosity_nu_start=luminosity_nu_start,
                   luminosity_nu_end=luminosity_nu_end,
                   last_no_of_packets=last_no_of_packets,
                   luminosity_requested=config.supernova.luminosity_requested.cgs,
                   convergence_strategy=config.montecarlo.convergence_strategy,
                   nthreads=config.montecarlo.nthreads)

