import warnings
import numpy as np
import pandas as pd
from astropy import units as u

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

        if not self.runner.line_interaction_type == 'downbranch':
            return raise_or_return(
                    'The FormalIntegrator currently only works for '
                    'line_interaction_type == "downbranch"'
                    )

        if self.runner.sigma_thomson.value > 1e-100:
            return raise_or_return(
                    'The FormalIntegrator currently only works with '
                    'disable_electron_scattering turned on.'
                    )

        return True

    def calculate_spectrum(self, frequency, points=None, raises=True):
        # Very crude implementation
        # The c extension needs bin centers (or something similar)
        # while TARDISSpectrum needs bin edges
        self.check(raises)
        N = points or self.points
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
        return result.ix[atomic_data.lines.index.values].as_matrix(), e_dot_u
