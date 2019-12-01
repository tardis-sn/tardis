import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.interpolate import interp1d
from astropy import units as u
from tardis import constants as const

from tardis.montecarlo.montecarlo import formal_integral
from tardis.montecarlo.spectrum import TARDISSpectrum


class IntegrationError(Exception):
    pass

class FormalIntegrator(object):

    def __init__(self, model, plasma, runner, points=1000):
        self.model = model
        self.plasma = plasma
        self.runner = runner
        self.points = points

    def check(self, raises=True):
        '''
        A method that determines if the formal integral can be performed with
        the current configuration settings

        The function returns False if the configuration conflicts with the
        required settings. If raises evaluates to True, then a
        IntegrationError is raised instead
        '''
        def raise_or_return(message):
            if raises:
                raise IntegrationError(message)
            else:
                warnings.warn(message)
                return False

        for obj in (self.model, self.plasma, self.runner):
            if obj is None:
                return raise_or_return(
                        'The integrator is missing either model, plasma or '
                        'runner. Please make sure these are provided to the '
                        'FormalIntegrator.'
                        )

        if not self.runner.line_interaction_type in ['downbranch', 'macroatom']:
            return raise_or_return(
                    'The FormalIntegrator currently only works for '
                    'line_interaction_type == "downbranch"'
                    'and line_interaction_type == "macroatom"'
                    )

        return True

    def calculate_spectrum(self, frequency, points=None,
            interpolate_shells=-1, raises=True):
        # Very crude implementation
        # The c extension needs bin centers (or something similar)
        # while TARDISSpectrum needs bin edges
        self.check(raises)
        N = points or self.points
        self.interpolate_shells = interpolate_shells
        frequency = frequency.to('Hz', u.spectral())

        luminosity = u.Quantity(
                formal_integral(
                    self,
                    frequency,
                    N),
                'erg'
                ) * (frequency[1] - frequency[0])

        # Ugly hack to convert to 'bin edges'
        frequency = u.Quantity(
                np.concatenate([
                    frequency.value,
                    [
                        frequency.value[-1] + np.diff(frequency.value)[-1]
                        ]]),
                    frequency.unit)

        return TARDISSpectrum(
                frequency,
                luminosity
                )

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
        macro_ref = atomic_data.macro_atom_references
        macro_data = atomic_data.macro_atom_data

        no_lvls = len(atomic_data.levels)
        no_shells = len(model.w)

        if runner.line_interaction_type == 'macroatom':
            internal_jump_mask = (macro_data.transition_type >= 0).values
            ma_int_data = macro_data[internal_jump_mask]
            internal = plasma.transition_probabilities[internal_jump_mask]

            source_level_idx = ma_int_data.source_level_idx.values
            destination_level_idx = ma_int_data.destination_level_idx.values 

        Edotlu_norm_factor = (1 / (runner.time_of_simulation * model.volume))
        exptau = 1 - np.exp(- plasma.tau_sobolevs)
        Edotlu = Edotlu_norm_factor * exptau * runner.Edotlu_estimator

        # The following may be achieved by calling the appropriate plasma
        # functions
        Jbluelu_norm_factor = (const.c.cgs * model.time_explosion /
                                (4 * np.pi * runner.time_of_simulation *
                                 model.volume)).to("1/(cm^2 s)").value
        # Jbluelu should already by in the correct order, i.e. by wavelength of
        # the transition l->u
        Jbluelu = runner.j_blue_estimator * Jbluelu_norm_factor

        upper_level_index = atomic_data.lines.index.droplevel('level_number_lower')
        e_dot_lu          = pd.DataFrame(Edotlu, index=upper_level_index)
        e_dot_u           = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u_src_idx = macro_ref.loc[e_dot_u.index].references_idx.values

        if runner.line_interaction_type == 'macroatom':
            C_frame = pd.DataFrame(
                columns=np.arange(no_shells), index=macro_ref.index
            )
            q_indices = (source_level_idx, destination_level_idx)
            for shell in range(no_shells):
                Q = sp.coo_matrix(
                        (internal[shell], q_indices), shape=(no_lvls, no_lvls)
                )
                inv_N = sp.identity(no_lvls) - Q
                e_dot_u_vec = np.zeros(no_lvls)
                e_dot_u_vec[e_dot_u_src_idx] = e_dot_u[shell].values
                C_frame[shell] = sp.linalg.spsolve(inv_N.T, e_dot_u_vec)

        e_dot_u.index.names = ['atomic_number', 'ion_number', 'source_level_number'] # To make the q_ul e_dot_u product work, could be cleaner
        transitions       = atomic_data.macro_atom_data[atomic_data.macro_atom_data.transition_type == -1].copy()
        transitions_index = transitions.set_index(['atomic_number', 'ion_number', 'source_level_number']).index.copy()
        tmp  = plasma.transition_probabilities[(atomic_data.macro_atom_data.transition_type == -1).values]
        q_ul = tmp.set_index(transitions_index)
        t    = model.time_explosion.value
        lines = atomic_data.lines.set_index('line_id')
        wave = lines.wavelength_cm.loc[transitions.transition_line_id].values.reshape(-1,1)
        if runner.line_interaction_type == 'macroatom':
            e_dot_u = C_frame.loc[e_dot_u.index]
        att_S_ul = (wave * (q_ul * e_dot_u) * t  / (4 * np.pi))

        result = pd.DataFrame(att_S_ul.values, index=transitions.transition_line_id.values)
        att_S_ul = result.loc[lines.index.values].values

        # Jredlu should already by in the correct order, i.e. by wavelength of
        # the transition l->u (similar to Jbluelu)
        Jredlu = Jbluelu * np.exp(-plasma.tau_sobolevs.values) + att_S_ul
        if self.interpolate_shells > 0:
            att_S_ul, Jredlu, Jbluelu, e_dot_u = self.interpolate_integrator_quantities(
                    att_S_ul, Jredlu, Jbluelu, e_dot_u)
        else:
            runner.r_inner_i = runner.r_inner_cgs
            runner.r_outer_i = runner.r_outer_cgs
            runner.tau_sobolevs_integ = plasma.tau_sobolevs.values
            runner.electron_densities_integ = plasma.electron_densities.values

        return att_S_ul, Jredlu, Jbluelu, e_dot_u

    def interpolate_integrator_quantities(self, att_S_ul, Jredlu,
            Jbluelu, e_dot_u):
        runner = self.runner
        plasma = self.plasma
        nshells = self.interpolate_shells
        r_middle = (runner.r_inner_cgs + runner.r_outer_cgs) / 2.

        r_integ = np.linspace(
                runner.r_inner_cgs[0], runner.r_outer_cgs[-1], nshells
        )
        runner.r_inner_i = r_integ[:-1]
        runner.r_outer_i = r_integ[1:]

        r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.

        runner.electron_densities_integ = interp1d(
                r_middle, plasma.electron_densities,
                fill_value='extrapolate', kind='nearest')(r_middle_integ)
        # Assume tau_sobolevs to be constant within a shell
        # (as in the MC simulation)
        runner.tau_sobolevs_integ = interp1d(
                r_middle, plasma.tau_sobolevs,
                fill_value='extrapolate', kind='nearest')(r_middle_integ)
        att_S_ul = interp1d(
                r_middle, att_S_ul, fill_value='extrapolate')(r_middle_integ)
        Jredlu = interp1d(
                r_middle, Jredlu, fill_value='extrapolate')(r_middle_integ)
        Jbluelu = interp1d(
                r_middle, Jbluelu, fill_value='extrapolate')(r_middle_integ)
        e_dot_u = interp1d(
                r_middle, e_dot_u, fill_value='extrapolate')(r_middle_integ)

        # Set negative values from the extrapolation to zero
        att_S_ul = att_S_ul.clip(0.)
        Jbluelu = Jbluelu.clip(0.)
        Jredlu = Jredlu.clip(0.)
        e_dot_u = e_dot_u.clip(0.)
        return att_S_ul, Jredlu, Jbluelu, e_dot_u
