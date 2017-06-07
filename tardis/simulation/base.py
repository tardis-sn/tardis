import time
import logging
import numpy as np
import pandas as pd
from astropy import units as u, constants as c
from collections import OrderedDict

from tardis.montecarlo import MontecarloRunner
from tardis.model import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.util import intensity_black_body

# Adding logging support
logger = logging.getLogger(__name__)


class Simulation(object):
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
    def __init__(self, iterations, model, plasma, runner,
                 no_of_packets, no_of_virtual_packets, luminosity_nu_start,
                 luminosity_nu_end, last_no_of_packets,
                 luminosity_requested, convergence_strategy,
                 nthreads):
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
        if convergence_strategy.type in ('damped', 'specific'):
            self.convergence_strategy = convergence_strategy
            self.converged = False
            self.consecutive_converges_count = 0
        else:
            raise ValueError('Convergence strategy type is '
                             'neither damped nor specific '
                             '- input is {0}'.format(convergence_strategy.type))

        self._callbacks = OrderedDict()
        self._cb_next_id = 0

    def estimate_t_inner(self, input_t_inner, luminosity_requested,
                         t_inner_update_exponent=-0.5):
        luminosity_ratios = (
            (self.emitted_luminosity / luminosity_requested).to(1).value)

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

        self.emitted_luminosity = self.runner.calculate_emitted_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end)
        reabsorbed_luminosity = self.runner.calculate_reabsorbed_luminosity(
            self.luminosity_nu_start, self.luminosity_nu_end)
        self.log_run_results(self.emitted_luminosity,
                             reabsorbed_luminosity)
        self.iterations_executed += 1

    def run(self):
        start_time = time.time()
        while self.iterations_executed < self.iterations-1 and not self.converged:
            self.iterate(self.no_of_packets)
            self.converged = self.advance_state()
            self._call_back()
        # Last iteration
        self.iterate(self.last_no_of_packets, self.no_of_virtual_packets, True)

        logger.info("Simulation finished in {0:d} iterations "
                    "and took {1:.2f} s".format(
                        self.iterations_executed, time.time() - start_time))
        self._call_back()
        # TODO: Add config flag to activate formal integral
        # self.runner.att_S_ul,self.runner.Edot_u =  self.make_source_function(model)
        # self.runner.L_nu,self.runner.L_nu_nus   =  self.integrate(model)

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

    def to_hdf(self, path_or_buf, path='simulation', plasma_properties=None,
               suffix_count=True):
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
        suffix_count : bool
            If True, the path inside the HDF will be suffixed with the
            number of the iteration being stored.
        Returns
        -------
        None
        """
        if suffix_count:
            path += str(self.iterations_executed)
        self.runner.to_hdf(path_or_buf, path)
        self.model.to_hdf(path_or_buf, path)
        self.plasma.to_hdf(path_or_buf, path, plasma_properties)

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

    def make_source_function(self):
        """
        Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

        Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.

        Parameters
        ----------
        model : tardis.model.Radial1DModel

        Returns
        -------
        Numpy array containing ( 1 - exp(-tau_ul) ) S_ul ordered by wavelength of the transition u -> l
        """
        model = self.model
        plasma = self.plasma
        runner = self.runner
        atomic_data = self.plasma.atomic_data

        Edotlu_norm_factor = (1 / (runner.time_of_simulation * model.volume))
        exptau = 1 - np.exp(- plasma.tau_sobolevs)
        Edotlu = Edotlu_norm_factor * exptau * runner.Edotlu_estimator

        upper_level_index = atomic_data.lines.set_index(['atomic_number', 'ion_number', 'level_number_upper']).index.copy()
        e_dot_lu          = pd.DataFrame(Edotlu, index=upper_level_index)
        e_dot_u           = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u.index.names = ['atomic_number', 'ion_number', 'source_level_number'] # To make the q_ul e_dot_u product work, could be cleaner
        transitions       = atomic_data.macro_atom_data[atomic_data.macro_atom_data.transition_type == -1].copy()
        transitions_index = transitions.set_index(['atomic_number', 'ion_number', 'source_level_number']).index.copy()
        tmp  = plasma.transition_probabilities[(atomic_data.macro_atom_data.transition_type == -1).values]
        q_ul = tmp.set_index(transitions_index)
        t    = model.time_explosion.value
        wave = atomic_data.lines.wavelength_cm[transitions.transition_line_id].values.reshape(-1,1)
        att_S_ul =  ( wave * (q_ul * e_dot_u) * t  / (4*np.pi) )

        result = pd.DataFrame(att_S_ul.as_matrix(), index=transitions.transition_line_id.values)
        return result.ix[atomic_data.lines.index.values].as_matrix(),e_dot_u

    def integrate(self):
        model = self.model
        plasma = self.plasma
        runner = self.runner

        num_shell, = self.runner.volume.shape
        ps         = np.linspace(0.999, 0, num = 20) # 3 * num_shell)
        R_max      = runner.r_outer_cgs.max()
        R_min_rel  = runner.r_inner_cgs.min() / R_max
        ct         = c.c.cgs.value * model.time_explosion.value / R_max

        att_S_ul, Edot_u = self.make_source_function()
        J_blues    = plasma.j_blues
        J_rlues    = plasma.j_blues * np.exp( -plasma.tau_sobolevs.values) + att_S_ul

        r_shells = np.zeros((num_shell+1,1))
        # Note the reorder from outer to inner
        r_shells[1:,0],r_shells[0,0] = self.runner.r_inner_cgs[::-1] / R_max, 1.0
        z_crossings = np.sqrt(r_shells**2 - ps**2)

        z_ct = z_crossings/ct
        z_ct[np.isnan(z_crossings)] = 0

        ## p > Rmin
        ps_outer        = ps[ ps > R_min_rel]
        z_ct_outer      = z_ct[:, ps > R_min_rel]
        n_shell_p_outer = (num_shell+1) - np.isnan(z_crossings[:, ps > R_min_rel]).sum(axis=0)

        ## p < Rmin
        ps_inner        = ps[ ps <= R_min_rel]
        z_ct_inner      = z_ct[:, ps <= R_min_rel]


        #Allocating stuff
        nus = runner.spectrum.frequency
        L_nu  = np.zeros(nus.shape)

        #Just aliasing for cleaner expressions later
        line_nu  = plasma.lines.nu
        taus     = plasma.tau_sobolevs
        T        = model.t_inner

        L_nu = self.py_integrate(L_nu, line_nu.as_matrix(), taus.as_matrix(), att_S_ul, R_max,
                    T.value, nus, ps_outer, ps_inner, z_ct_outer, z_ct_inner,
                    num_shell,n_shell_p_outer)

        return  L_nu#,nus

    def py_integrate(self, L_nu, line_nu, taus, att_S_ul, R_max,
                    T, nus, ps_outer, ps_inner, z_ct_outer, z_ct_inner,
                    num_shell,n_shell_p_outer):
    # def py_integrate(self,L_nu, line_nu, taus, att_S_ul, T, nus, num_shell, R_max, ps_outer, ps_inner, n_shell_p_outer, z_ct_outer, z_ct_inner):
        for nu_idx,nu in enumerate(nus.value):
            I_outer = np.zeros(ps_outer.shape)
            for p_idx,p in enumerate(ps_outer):
                z_cross_p = z_ct_outer[z_ct_outer[:,p_idx] > 0,p_idx]
                if len(z_cross_p) == 1:
                    z_cross_p = np.hstack((-z_cross_p,z_cross_p))
                else:
                    z_cross_p = np.hstack((-z_cross_p,z_cross_p[::-1][1:],0)) # Zero ensures empty ks in last step below
                                                                              # 1: avoids double counting center
                shell_idx = (num_shell-1) - np.arange(n_shell_p_outer[p_idx]) # -1 for 0-based indexing
                shell_idx = np.hstack((shell_idx,shell_idx[::-1][1:]))

                for idx,z_cross in enumerate(z_cross_p[:-1]):
                    if z_cross_p[idx+1] == 0:
                        continue
                    nu_start = nu * (1 - z_cross)
                    nu_end   = nu * (1 - z_cross_p[idx+1])
                    shell = shell_idx[idx]
                    # Note the direction of the comparisons
                    ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) )

                    if len(ks) < 2:
                        ks = list(ks)

                    #I_outer[p_idx] = I_outer[p_idx] + dtau * (
                    #                    ( J_rlues.iloc[ks[0],shell] + J_blues.iloc[ks[1],shell] ) / 2 - I_outer[p_idx] )
                    for k in ks:
                        I_outer[p_idx] = I_outer[p_idx] * np.exp(-taus[k,shell]) + att_S_ul[k,shell]


            I_inner = np.zeros(ps_inner.shape)
            for p_idx,p in enumerate(ps_inner):
                z_cross_p = z_ct_inner[z_ct_inner[:,p_idx] > 0,p_idx]
                z_cross_p = np.hstack((z_cross_p[::-1],0)) # Zero ensures empty ks in last step below

                shell_idx = np.hstack(( np.arange(num_shell), 0 ))
                I_inner[p_idx] = intensity_black_body(nu,T)
                for idx,z_cross in enumerate(z_cross_p[:-1]):
                    if z_cross_p[idx+1] == 0:
                        continue
                    nu_start = nu * (1 - z_cross)
                    nu_end   = nu * (1 - z_cross_p[idx+1])
                    shell = shell_idx[idx]
                    # Note the direction of the comparisons
                    ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) )

                    if len(ks) < 2:
                        ks = list(ks)

                    #dtau * (
                    #( J_rlues.iloc[ks[0],shell] + J_blues.iloc[ks[1],shell] ) / 2 )

                    for k in ks:
                        I_inner[p_idx] = I_inner[p_idx] * np.exp(-taus[k,shell]) + att_S_ul[k,shell]
                        if ( nu_idx % 200 ) == 0:
                            print p_idx,idx, k, I_inner[p_idx]
                    if (( nu_idx % 200 ) == 0) & (idx == 1):
                        print p_idx,idx, I_inner[p_idx]

            if ( nu_idx % 200 ) == 0:
                print "{:3.0f} %".format( 100*float(nu_idx)/len(nus))
                print I_outer, I_inner
            ps = np.hstack((ps_outer,ps_inner))*R_max
            I_p = np.hstack((I_outer,I_inner))*ps
            L_nu[nu_idx] = 8 * np.pi**2 *  np.trapz(y = I_p[::-1],x = ps[::-1])

        return L_nu
