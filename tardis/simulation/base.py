import os
import logging
import time
import itertools

import pandas as pd

import numpy as np

from astropy import units as u
from astropy import constants as co

from tardis.montecarlo import montecarlo
from tardis.montecarlo.base import MontecarloRunner
from tardis.plasma.properties.base import Input
from tardis.util import intensity_black_body
# Adding logging support
logger = logging.getLogger(__name__)


class Simulation(object):

    converged = False

    def __init__(self, tardis_config):
        self.tardis_config = tardis_config
        self.runner = MontecarloRunner(self.tardis_config.montecarlo.seed,
                                       tardis_config.spectrum.frequency,
                                       tardis_config.supernova.get('distance',
                                                                   None))
        t_inner_lock_cycle = [False] * (tardis_config.montecarlo.
                                        convergence_strategy.
                                        lock_t_inner_cycles)
        t_inner_lock_cycle[0] = True
        self.t_inner_update = itertools.cycle(t_inner_lock_cycle)

    def run_single_montecarlo(self, model, no_of_packets,
                              no_of_virtual_packets=0,last_run=False):
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
                        nthreads=self.tardis_config.montecarlo.nthreads,last_run=last_run)


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

    def legacy_run_simulation(self, model, hdf_path_or_buf=None,
                              hdf_mode='full', hdf_last_only=True):
        """

        Parameters
        ----------
        model : tardis.model.Radial1DModel
        hdf_path_or_buf : str, optional
            A path to store the data of each simulation iteration
            (the default value is None, which means that nothing
            will be stored).
        hdf_mode : {'full', 'input'}, optional
            If 'full' all plasma properties will be stored to HDF,
            if 'input' only input plasma properties will be stored.
        hdf_last_only: bool, optional
            If True, only the last iteration of the simulation will
            be stored to the HDFStore.

        Returns
        -------

        """
        if hdf_path_or_buf is not None:
            if hdf_mode == 'full':
                plasma_properties = None
            elif hdf_mode == 'input':
                plasma_properties = [Input]
            else:
                raise ValueError('hdf_mode must be "full" or "input"'
                                 ', not "{}"'.format(type(hdf_mode)))
        start_time = time.time()

        self.iterations_remaining = self.tardis_config.montecarlo.iterations
        self.iterations_max_requested = self.tardis_config.montecarlo.iterations
        self.iterations_executed = 0
        converged = False

        convergence_section = (
                    self.tardis_config.montecarlo.convergence_strategy)

        while self.iterations_remaining > 1:
            logger.info('Remaining run %d', self.iterations_remaining)
            self.run_single_montecarlo(
                model, self.tardis_config.montecarlo.no_of_packets)
            self.log_run_results(self.calculate_emitted_luminosity(),
                                 self.calculate_reabsorbed_luminosity())
            self.iterations_executed += 1
            self.iterations_remaining -= 1

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
            if hdf_path_or_buf is not None and not hdf_last_only:
                self.to_hdf(model, hdf_path_or_buf,
                            'simulation{}'.format(self.iterations_executed),
                            plasma_properties)


            # if switching into the hold iterations mode or out back to the normal one
            # if it is in either of these modes already it will just stay there
            if converged and not self.converged:
                self.converged = True
                # UMN - used to be 'hold_iterations_wrong' but this is
                # currently not in the convergence_section namespace...
                self.iterations_remaining = (
                    convergence_section["hold_iterations"])
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
                self.iterations_remaining = (
                    convergence_section["hold_iterations"])

        #Finished second to last loop running one more time
        logger.info('Doing last run')
        if self.tardis_config.montecarlo.last_no_of_packets is not None:
            no_of_packets = self.tardis_config.montecarlo.last_no_of_packets
        else:
            no_of_packets = self.tardis_config.montecarlo.no_of_packets

        no_of_virtual_packets = (
            self.tardis_config.montecarlo.no_of_virtual_packets)

        self.run_single_montecarlo(model, no_of_packets, no_of_virtual_packets, last_run=True)

        self.runner.legacy_update_spectrum(no_of_virtual_packets)
        self.legacy_set_final_model_properties(model)

        #the following instructions, passing down information to the model are
        #required for the gui
        model.no_of_packets = no_of_packets
        model.no_of_virtual_packets = no_of_virtual_packets
        model.converged = converged
        model.iterations_executed = self.iterations_executed
        model.iterations_max_requested = self.iterations_max_requested

        logger.info("Finished in {0:d} iterations and took {1:.2f} s".format(
            self.iterations_executed, time.time()-start_time))

        if hdf_path_or_buf is not None:
            self.to_hdf(model, hdf_path_or_buf,
                        'simulation{}'.format(self.iterations_executed),
                        plasma_properties)

        self.runner.att_S_ul,self.runner.Edot_u,self.runner.wave =  self.make_source_function(model)
        self.runner.L_nu,self.runner.L_nu_nus     =  self.integrate(model)

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

    def to_hdf(self, model, path_or_buf, path='', plasma_properties=None):
        """
        Store the simulation to an HDF structure.

        Parameters
        ----------
        model : tardis.model.Radial1DModel
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
        model.to_hdf(path_or_buf, path, plasma_properties)

    def make_source_function(self, model):
        """
        Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

        Formally it calculates the expresion ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.

        Parameters
        ----------
        model : tardis.model.Radial1DModel

        Returns
        -------
        DataFrame containing ( 1 - exp(-tau_ul) ) S_ul
        """

        upper_level_index = model.atom_data.lines.set_index(['atomic_number', 'ion_number', 'level_number_upper']).index.copy()
        e_dot_lu          = pd.DataFrame(self.runner.Edotlu, index=upper_level_index)
        e_dot_u           = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u.index.names = ['atomic_number', 'ion_number', 'source_level_number'] # To make the q_ul e_dot_u product work, could be cleaner
        transitions       = model.atom_data.macro_atom_data[model.atom_data.macro_atom_data.transition_type == -1].copy()
        transitions_index = transitions.set_index(['atomic_number', 'ion_number', 'source_level_number']).index.copy()
        tmp  = model.plasma.transition_probabilities[(model.atom_data.macro_atom_data.transition_type == -1).values]
        q_ul = tmp.set_index(transitions_index)
        t    = model.tardis_config.supernova.time_explosion.value
        wave = model.atom_data.lines.wavelength_cm[transitions.transition_line_id].values.reshape(-1,1)
        reorder  = wave[:,0].argsort()
        att_S_ul =  ( wave * (q_ul * e_dot_u) * t  / (4*np.pi) ).iloc[reorder,:]

        return att_S_ul.as_matrix(),e_dot_u,wave

    def integrate(self,model):
        num_shell, = self.runner.volume.shape
        ps         = np.linspace(0.999, 0, num = 10) # 3 * num_shell)
        R_max      = self.runner.r_outer_cgs.max()
        R_min_rel  = self.runner.r_inner_cgs.min() / R_max
        ct         = co.c.cgs.value * self.tardis_config.supernova.time_explosion.value / R_max
        J_blues    = model.j_blues
#        J_rlues    = model.j_blues.values * np.exp( -model.plasma.tau_sobolevs.values) + self.runner.att_S_ul

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

        # I will traverse from largest shell in, but
        # elsewhere shell structure is from smallest and out,
        # so this is for reversing it
        shell_nr = np.arange(0,num_shell,dtype="int")[::-1]

 
        nus = self.runner.spectrum.frequency


        #Just aliasing for cleaner expressions later
        line_nu  = model.plasma.lines.nu
        taus     = model.plasma.tau_sobolevs
        att_S_ul = self.runner.att_S_ul
        T        = model.t_inner

        dtau = 0 # Just to remember it 
        cnst = 0

        L_nu  = np.zeros(nus.shape)

        #### Debug ####
        import pdb 

        for nu_idx,nu in enumerate(nus.value):
            I_outer = np.zeros(ps_outer.shape)
            for p_idx,p in enumerate(ps_outer):
                z_cross_p = z_ct_outer[z_ct_outer[:,p_idx] > 0,p_idx]
                z_cross_p = np.hstack((-z_cross_p,z_cross_p[::-1][1:],0)) # Zero ensures empty ks in last step below
                                                                          # 1: avoids double counting center
                shell_idx = (num_shell-1) - np.arange(n_shell_p_outer[p_idx]) # -1 for 0-based indexing
                shell_idx = np.hstack((shell_idx,shell_idx[::-1][1:]))
                
                for idx,z_cross in enumerate(z_cross_p[:-1]):
                    nu_start = nu / (1 + z_cross) 
                    nu_end   = nu / (1 + z_cross_p[idx+1])
                    shell = shell_idx[idx]
                    # Note the direction of the comparisons
                    ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) ) 

                    if len(ks) < 2:
                        continue

                    #I_outer[p_idx] = I_outer[p_idx] + dtau * ( 
                    #                    ( J_rlues.iloc[ks[0],shell] + J_blues.iloc[ks[1],shell] ) / 2 - I_outer[p_idx] )
                    for k in ks:
                        I_outer[p_idx] = I_outer[p_idx] * np.exp(-taus.iloc[k,shell]) + att_S_ul[k,shell]


            I_inner = np.zeros(ps_inner.shape)
            for p_idx,p in enumerate(ps_inner):
                z_cross_p = z_ct_inner[z_ct_inner[:,p_idx] > 0,p_idx]
                z_cross_p = np.hstack((z_cross_p[::-1],0)) # Zero ensures empty ks in last step below

                shell_idx = np.hstack(( np.arange(num_shell), 0 ))
                I_inner[p_idx] = intensity_black_body(nu,T)
                for idx,z_cross in enumerate(z_cross_p[:-1]):
                    nu_start = nu / (1 + z_cross) 
                    nu_end   = nu / (1 + z_cross_p[idx+1])
                    shell = shell_idx[idx]
                    # Note the direction of the comparisons
                    ks, = np.where( (line_nu < nu_start) & (line_nu >= nu_end) )                    

                    if len(ks) < 2:
                        continue


                    #dtau * ( 
                    #( J_rlues.iloc[ks[0],shell] + J_blues.iloc[ks[1],shell] ) / 2 )

#                    if p_idx == 0:
#                        print "p_idx", p_idx

                    for k in ks:
                        I_inner[p_idx] = I_inner[p_idx] * np.exp(-taus.iloc[k,shell]) + att_S_ul[k,shell]
#                        if (p_idx == 0) &  (I_inner[p_idx] > 0.0001):
#                            print  att_S_ul.iloc[k,shell],np.exp(-taus.iloc[k,shell])


            if ( nu_idx % 10 ) == 0:
                print "{:3.0f} %".format( 100*float(nu_idx)/len(nus))
                print I_outer, I_inner
            ps = np.hstack((ps_outer,ps_inner))*R_max
            I_p = np.hstack((I_outer,I_inner))*ps
            L_nu[nu_idx] = 8 * np.pi**2 *  np.trapz(y = I_p[::-1],x = ps[::-1])

        return  L_nu,nus


def run_radial1d(radial1d_model, hdf_path_or_buf=None,
                 hdf_mode='full', hdf_last_only=True):
    """

    Parameters
    ----------
    radial1d_model : tardis.model.Radial1DModel
    hdf_path_or_buf : str, optional
        A path to store the data of each simulation iteration
        (the default value is None, which means that nothing
        will be stored).
    hdf_mode : {'full', 'input'}, optional
        If 'full' all plasma properties will be stored to HDF,
        if 'input' only input plasma properties will be stored.
    hdf_last_only: bool, optional
        If True, only the last iteration of the simulation will
        be stored to the HDFStore.

    Returns
    -------

    """

    simulation = Simulation(radial1d_model.tardis_config)
    simulation.legacy_run_simulation(radial1d_model, hdf_path_or_buf,
                                     hdf_mode, hdf_last_only)
