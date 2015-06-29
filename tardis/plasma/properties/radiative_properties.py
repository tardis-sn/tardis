import logging

import numpy as np
import pandas as pd
from astropy import units as u, constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis import macro_atom

logger = logging.getLogger(__name__)

__all__ = ['StimulatedEmissionFactor', 'TauSobolev', 'BetaSobolev',
    'TransitionProbabilities']

class StimulatedEmissionFactor(ProcessingPlasmaProperty):

    outputs = ('stimulated_emission_factor',)

    def __init__(self, plasma_parent):
        super(StimulatedEmissionFactor, self).__init__(plasma_parent)
        self._g_upper = None
        self._g_lower = None

    def get_g_lower(self, levels, lines_lower_level_index):
        if self._g_lower is None:
            g_lower = np.array(levels.g.ix[lines_lower_level_index],
                                     dtype=np.float64)
            self._g_lower = g_lower[np.newaxis].T
        return self._g_lower


    def get_g_upper(self, levels, lines_upper_level_index):
        if self._g_upper is None:
            g_upper = np.array(levels.g.ix[lines_upper_level_index],
                                     dtype=np.float64)
            self._g_upper = g_upper[np.newaxis].T
        return self._g_upper

    def calculate(self, levels, level_number_density, lines_lower_level_index,
        lines_upper_level_index):

        n_lower = level_number_density.values.take(lines_lower_level_index,
            axis=0, mode='raise').copy('F')
        n_upper = level_number_density.values.take(lines_upper_level_index,
            axis=0, mode='raise').copy('F')

        meta_stable_upper = levels.metastable.values.take(
            lines_upper_level_index, axis=0, mode='raise')[np.newaxis].T

        g_lower = self.get_g_lower(levels, lines_lower_level_index)
        g_upper = self.get_g_upper(levels, lines_upper_level_index)

        stimulated_emission_factor = 1 - ((g_lower * n_upper) / (g_upper * n_lower))

        stimulated_emission_factor[n_lower == 0.0] = 0.0
        stimulated_emission_factor[np.isneginf(stimulated_emission_factor)] = 0.0
        stimulated_emission_factor[meta_stable_upper &
                                   (stimulated_emission_factor < 0)] = 0.0

        return stimulated_emission_factor


class TauSobolev(ProcessingPlasmaProperty):
    """
    This function calculates the Sobolev optical depth :math:`\\tau_\\textrm{Sobolev}`

    .. math::
        C_\\textrm{Sobolev} = \\frac{\\pi e^2}{m_e c}

        \\tau_\\textrm{Sobolev} = C_\\textrm{Sobolev}\,  \\lambda\\, f_{\\textrm{lower}\\rightarrow\\textrm{upper}}\\,
            t_\\textrm{explosion}\, N_\\textrm{lower}



    .. note::
        Currently we're ignoring the term for stimulated emission:
            :math:`(1 - \\frac{g_\\textrm{lower}}{g_\\textrm{upper}}\\frac{N_\\textrm{upper}}{N_\\textrm{lower}})`


    """

    outputs = ('tau_sobolevs',)


    def __init__(self, plasma_parent):
        super(TauSobolev, self).__init__(plasma_parent)
        self.sobolev_coefficient = (((np.pi * const.e.gauss ** 2) /
                                    (const.m_e.cgs * const.c.cgs))
                                    * u.cm * u.s / u.cm**3).to(1).value

    def calculate(self, lines, level_number_density, lines_lower_level_index,
                  time_explosion, stimulated_emission_factor, j_blues):

        f_lu = lines.f_lu.values[np.newaxis].T
        wavelength = lines.wavelength_cm.values[np.newaxis].T

        n_lower = level_number_density.values.take(lines_lower_level_index, axis=0, mode='raise').copy('F')

        #if self.nlte_config is not None and self.nlte_config.species != []:
        #    nlte_lines_mask = np.zeros(self.stimulated_emission_factor.shape[0]).astype(bool)
        #    for species in self.nlte_config.species:
        #        nlte_lines_mask |= (self.atom_data.lines.atomic_number == species[0]) & \
        #                           (self.atom_data.lines.ion_number == species[1])
        #    self.stimulated_emission_factor[(self.stimulated_emission_factor < 0) & nlte_lines_mask[np.newaxis].T] = 0.0


        tau_sobolevs = (self.sobolev_coefficient * f_lu * wavelength *
                        time_explosion * n_lower * stimulated_emission_factor)

        return pd.DataFrame(tau_sobolevs, index=lines.index,
            columns=np.array(level_number_density.columns))

class BetaSobolev(ProcessingPlasmaProperty):
    outputs = ('beta_sobolev',)

    def calculate(self, tau_sobolevs):
        if not hasattr(self, 'beta_sobolev'):
            beta_sobolev = np.zeros_like(tau_sobolevs.values)
        macro_atom.calculate_beta_sobolev(
            tau_sobolevs.values.ravel(order='F'),
            beta_sobolev.ravel(order='F'))
        return beta_sobolev

class TransitionProbabilities(ProcessingPlasmaProperty):
    outputs = ('transition_probabilities',)

    def calculate(self, atomic_data, beta_sobolev, j_blues,
        stimulated_emission_factor, tau_sobolevs):
        if j_blues.empty:
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
            j_blues = j_blues.values.take(macro_atom_transition_up_filter,
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
