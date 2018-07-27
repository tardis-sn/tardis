import time
import logging
import numpy as np
import pandas as pd
from astropy import units as u
from collections import OrderedDict

from tardis.montecarlo import MontecarloRunner
from tardis.model import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.io.util import HDFWriterMixin
# Adding logging support
logger = logging.getLogger(__name__)


class PlasmaStateStorerMixin(object):
    """Mixin class to provide the capability to the simulation object of
    storing plasma information and the inner boundary temperature during each
    MC iteration.

    Currently, storage for the dilution factor, the radiation temperature and
    the electron density in each cell is provided. Additionally, the
    temperature at the inner boundary is saved.
    """
    def __init__(self, iterations, no_of_shells):

        self.iterations_w = np.zeros(
            (iterations, no_of_shells))
        self.iterations_t_rad = np.zeros(
            (iterations, no_of_shells)) * u.K
        self.iterations_electron_densities = np.zeros(
            (iterations, no_of_shells))
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
        self.iterations_electron_densities[i, :] = \
            electron_densities.values
        self.iterations_t_inner[i] = t_inner

    def reshape_plasma_state_store(self, executed_iterations):
        """Reshapes the storage arrays in case convergence was reached before
        all specified iterations were executed.

        Parameters
        ----------
        executed_iterations : int
            iteration index, i.e. number of iterations executed minus one!
        """
        self.iterations_w = self.iterations_w[:executed_iterations+1, :]
        self.iterations_t_rad = \
            self.iterations_t_rad[:executed_iterations+1, :]
        self.iterations_electron_densities = \
            self.iterations_electron_densities[:executed_iterations+1, :]
        self.iterations_t_inner = \
            self.iterations_t_inner[:executed_iterations+1]


class Simulation(PlasmaStateStorerMixin, HDFWriterMixin):
    """A composite object containing all the required information for a
    simulation.

    Parameters
    ----------
    converged : bool
    iterations : int
    model : tardis.model.Radial1DModel
    plasma : tardis.plasma.BasePlasma
    runner : tardis.montecarlo.MontecarloRunner
    no_of_packets : int
    last_no_of_packets : int
    no_of_virtual_packets : int
    luminosity_nu_start : astropy.units.Quantity
    luminosity_nu_end : astropy.units.Quantity
    luminosity_requested : astropy.units.Quantity
    nthreads : int
        The number of threads to run montecarlo with

        .. note:: TARDIS must be built with OpenMP support in order for
        `nthreads` to have effect.
    """
    hdf_properties = ['model', 'plasma', 'runner', 'iterations_w',
                      'iterations_t_rad', 'iterations_electron_densities',
                      'iterations_t_inner']
    hdf_name = 'simulation'
    def __init__(self, iterations, model, plasma, runner,
                 no_of_packets, no_of_virtual_packets, luminosity_nu_start,
                 luminosity_nu_end, last_no_of_packets,
                 luminosity_requested, convergence_strategy,
                 nthreads):

        super(Simulation, self).__init__(iterations, model.no_of_shells)

        self.converged = False
        self.iterations = iterations
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
        if convergence_strategy.type in ('damped'):
            self.convergence_strategy = convergence_strategy
            self.converged = False
            self.consecutive_converges_count = 0
        elif convergence_strategy.type in ('custom'):
            raise NotImplementedError(
                'Convergence strategy type is custom; '
                'you need to implement your specific treatment!'
            )
        else:
            raise ValueError(
                    'Convergence strategy type is '
                    'not damped or custom '
                    '- input is {0}'.format(convergence_strategy.type))

        self._callbacks = OrderedDict()
        self._cb_next_id = 0

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

        fraction_t_rad_converged = (
            np.count_nonzero(
                convergence_t_rad < self.convergence_strategy.t_rad.threshold)
            / no_of_shells)

        t_rad_converged = (
            fraction_t_rad_converged > self.convergence_strategy.fraction)

        fraction_w_converged = (
            np.count_nonzero(
                convergence_w < self.convergence_strategy.w.threshold)
            / no_of_shells)

        w_converged = (
            fraction_w_converged > self.convergence_strategy.fraction)

        t_inner_converged = (
            convergence_t_inner < self.convergence_strategy.t_inner.threshold)

        if np.all([t_rad_converged, w_converged, t_inner_converged]):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info("Iteration converged {0:d}/{1:d} consecutive "
                        "times.".format(self.consecutive_converges_count,
                                        hold_iterations + 1))
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
            converged : ~bool
        """
        estimated_t_rad, estimated_w = (
            self.runner.calculate_radiationfield_properties())
        estimated_t_inner = self.estimate_t_inner(
            self.model.t_inner, self.luminosity_requested,
            t_inner_update_exponent=self.convergence_strategy.t_inner_update_exponent)

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
        if (self.iterations_executed + 1) % self.convergence_strategy.lock_t_inner_cycles == 0:
            next_t_inner = self.damped_converge(
                self.model.t_inner, estimated_t_inner,
                self.convergence_strategy.t_inner.damping_constant)
        else:
            next_t_inner = self.model.t_inner

        self.log_plasma_state(self.model.t_rad, self.model.w,
                              self.model.t_inner, next_t_rad, next_w,
                              next_t_inner)
        self.model.t_rad = next_t_rad
        self.model.w = next_w
        self.model.t_inner = next_t_inner

        # model.calculate_j_blues() equivalent
        # model.update_plasmas() equivalent
        # Bad test to see if this is a nlte run
        if 'nlte_data' in self.plasma.outputs_dict:
            self.plasma.store_previous_properties()

        update_properties = dict(t_rad=self.model.t_rad, w=self.model.w)
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.
        if 'j_blue_estimator' in self.plasma.outputs_dict:
            update_properties.update(t_inner=next_t_inner,
                                 j_blue_estimator=self.runner.j_blue_estimator)

        self.plasma.update(**update_properties)

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
        while self.iterations_executed < self.iterations-1:
            self.store_plasma_state(self.iterations_executed, self.model.w,
                                    self.model.t_rad,
                                    self.plasma.electron_densities,
                                    self.model.t_inner)
            self.iterate(self.no_of_packets)
            self.converged = self.advance_state()
            self._call_back()
            if self.converged:
                if self.convergence_strategy.stop_if_converged:
                    break
        # Last iteration
        self.store_plasma_state(self.iterations_executed, self.model.w,
                                self.model.t_rad,
                                self.plasma.electron_densities,
                                self.model.t_inner)
        self.iterate(self.last_no_of_packets, self.no_of_virtual_packets, True)

        self.reshape_plasma_state_store(self.iterations_executed)

        logger.info("Simulation finished in {0:d} iterations "
                    "and took {1:.2f} s".format(
                        self.iterations_executed, time.time() - start_time))
        self._call_back()


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
        cb_func: callable
            The callback function
        arg1:
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
        id: int
            The callback ID

        Returns
        -------
        : True if the callback was successfully removed.
        """
        try:
            del self._callbacks[id]
            return True
        except KeyError:
            return False

    @classmethod
    def from_config(cls, config, **kwargs):
        """
        Create a new Simulation instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        **kwargs
            Allow overriding some structures, such as model, plasma, atomic data
            and the runner, instead of creating them from the configuration
            object.

        Returns
        -------
        Simulation

        """
        # Allow overriding some config structures. This is useful in some
        # unit tests, and could be extended in all the from_config classmethods.
        if 'model' in kwargs:
            model = kwargs['model']
        else:
            model = Radial1DModel.from_config(config)
        if 'plasma' in kwargs:
            plasma = kwargs['plasma']
        else:
            plasma = assemble_plasma(config, model,
                                     atom_data=kwargs.get('atom_data', None))
        if 'runner' in kwargs:
            runner = kwargs['runner']
        else:
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
