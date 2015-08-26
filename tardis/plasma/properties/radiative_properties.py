import logging

import numpy as np
import pandas as pd
from astropy import units as u, constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis import macro_atom

logger = logging.getLogger(__name__)

__all__ = ['StimulatedEmissionFactor', 'TauSobolev', 'BetaSobolev',
    'TransitionProbabilities', 'LTEJBlues']

class StimulatedEmissionFactor(ProcessingPlasmaProperty):
    """
    Outputs:
        stimulated_emission_factor : Numpy Array [len(lines), len(t_rad)]
    """
    outputs = ('stimulated_emission_factor',)
    latex_formula = ('1-\\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}',)

    def __init__(self, plasma_parent):
        super(StimulatedEmissionFactor, self).__init__(plasma_parent)
        self._g_upper = None
        self._g_lower = None

    def get_g_lower(self, g, lines_lower_level_index):
        if self._g_lower is None:
            g_lower = np.array(g.ix[lines_lower_level_index],
                                     dtype=np.float64)
            self._g_lower = g_lower[np.newaxis].T
        return self._g_lower

    def get_g_upper(self, g, lines_upper_level_index):
        if self._g_upper is None:
            g_upper = np.array(g.ix[lines_upper_level_index],
                                     dtype=np.float64)
            self._g_upper = g_upper[np.newaxis].T
        return self._g_upper

    def calculate(self, g, level_number_density, lines_lower_level_index,
        lines_upper_level_index, metastability, nlte_excitation_species, lines):
        n_lower = level_number_density.values.take(lines_lower_level_index,
            axis=0, mode='raise').copy('F')
        n_upper = level_number_density.values.take(lines_upper_level_index,
            axis=0, mode='raise').copy('F')
        meta_stable_upper = metastability.values.take(
            lines_upper_level_index, axis=0, mode='raise')[np.newaxis].T
        g_lower = self.get_g_lower(g, lines_lower_level_index)
        g_upper = self.get_g_upper(g, lines_upper_level_index)
        stimulated_emission_factor = 1 - ((g_lower * n_upper) /
            (g_upper * n_lower))
        stimulated_emission_factor[n_lower == 0.0] = 0.0
        stimulated_emission_factor[np.isneginf(stimulated_emission_factor)]\
            = 0.0
        stimulated_emission_factor[meta_stable_upper &
                                   (stimulated_emission_factor < 0)] = 0.0
        if nlte_excitation_species:
            nlte_lines_mask = \
                np.zeros(stimulated_emission_factor.shape[0]).astype(bool)
            for species in nlte_excitation_species:
                nlte_lines_mask |= (lines.atomic_number == species[0]) & \
                                   (lines.ion_number == species[1])
            stimulated_emission_factor[(stimulated_emission_factor < 0) &
                nlte_lines_mask[np.newaxis].T] = 0.0
        return stimulated_emission_factor

class TauSobolev(ProcessingPlasmaProperty):
    """
    Outputs:
        tau_sobolev : Pandas DataFrame
        Sobolev optical depth for each line.
    """
    outputs = ('tau_sobolevs',)
    latex_name = ('\\tau_{\\textrm{sobolev}}',)
    latex_formula = ('\\dfrac{\\pi e^{2}}{m_{e} c}f_{lu}\\lambda t_{exp}\
        n_{lower} \\Big(1-\\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\\Big)',)

    def __init__(self, plasma_parent):
        super(TauSobolev, self).__init__(plasma_parent)
        self.sobolev_coefficient = (((np.pi * const.e.gauss ** 2) /
                                    (const.m_e.cgs * const.c.cgs))
                                    * u.cm * u.s / u.cm**3).to(1).value

    def calculate(self, lines, level_number_density, lines_lower_level_index,
                  time_explosion, stimulated_emission_factor, j_blues,
                  f_lu, wavelength_cm):
        f_lu = f_lu.values[np.newaxis].T
        wavelength = wavelength_cm.values[np.newaxis].T
        n_lower = level_number_density.values.take(lines_lower_level_index,
            axis=0, mode='raise').copy('F')
        tau_sobolevs = (self.sobolev_coefficient * f_lu * wavelength *
                        time_explosion * n_lower * stimulated_emission_factor)
        return pd.DataFrame(tau_sobolevs, index=lines.index,
            columns=np.array(level_number_density.columns))

class BetaSobolev(ProcessingPlasmaProperty):
    """
    Outputs:
        beta_sobolev : Numpy Array
    """
    outputs = ('beta_sobolev',)
    latex_name = ('\\beta_{\\textrm{sobolev}}',)

    def calculate(self, tau_sobolevs):
        if getattr(self, 'beta_sobolev', None) is None:
            beta_sobolev = np.zeros_like(tau_sobolevs.values)
        else:
            beta_sobolev = self.beta_sobolev
        macro_atom.calculate_beta_sobolev(
            tau_sobolevs.values.ravel(order='F'),
            beta_sobolev.ravel(order='F'))
        return beta_sobolev

class TransitionProbabilities(ProcessingPlasmaProperty):
    """
    Outputs:
        transition_probabilities : Pandas DataFrame
    """
    outputs = ('transition_probabilities',)

    def calculate(self, atomic_data, beta_sobolev, j_blues,
        stimulated_emission_factor, tau_sobolevs):
        if len(j_blues) == 0:
            transition_probabilities = None
        else:
            try:
                macro_atom_data = atomic_data.macro_atom_data
            except:
                macro_atom_data = atomic_data.macro_atom_data_all
            transition_probabilities = (
                macro_atom_data.transition_probability.values[np.newaxis].T *
                beta_sobolev.take(macro_atom_data.lines_idx.values.astype(int),
                    axis=0, mode='raise')).copy('F')
            transition_up_filter = \
                (macro_atom_data.transition_type == 1).values
            macro_atom_transition_up_filter = \
                macro_atom_data.lines_idx.values[transition_up_filter]
            j_blues = j_blues.take(macro_atom_transition_up_filter,
                axis=0, mode='raise')
            macro_stimulated_emission = stimulated_emission_factor.take(
                macro_atom_transition_up_filter, axis=0, mode='raise')
            transition_probabilities[transition_up_filter] *= j_blues * \
                macro_stimulated_emission
            block_references = np.hstack((
                atomic_data.macro_atom_references.block_references,
                len(macro_atom_data)))
            macro_atom.normalize_transition_probabilities(
                transition_probabilities, block_references)
            transition_probabilities = pd.DataFrame(transition_probabilities,
                index=macro_atom_data.transition_line_id,
                columns=tau_sobolevs.columns)
        return transition_probabilities

class LTEJBlues(ProcessingPlasmaProperty):
    outputs = ('lte_j_blues',)
    latex_name = ('J^{b}_{lu(LTE)}')

    @staticmethod
    def calculate(lines, nu, beta_rad):
        beta_rad = pd.Series(beta_rad)
        nu = pd.Series(nu)
        h = const.h.cgs.value
        c = const.c.cgs.value
        df = pd.DataFrame(1, index=nu.index, columns=beta_rad.index)
        df = df.mul(nu, axis='index') * beta_rad
        exponential = (np.exp(h * df) - 1)**(-1)
        remainder = (2 * (h * nu.values ** 3) /
            (c ** 2))
        j_blues = exponential.mul(remainder, axis=0)
        return pd.DataFrame(j_blues, index=lines.index, columns=beta_rad.index)
