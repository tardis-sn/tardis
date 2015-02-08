import logging

import numpy as np
import pandas as pd

from tardis.plasma.plasma_properties import BasePlasmaProperty

logger = logging.getLogger(__name__)


class TauSobolev(BasePlasmaProperty):
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

    def calculate(self, lines):
        f_lu = self.atom_data.lines['f_lu'].values
        wavelength = self.atom_data.lines['wavelength_cm'].values


        #todo fix this is a concern the mode='safe'
        n_lower = self.level_populations.values.take(self.atom_data.lines_lower2level_idx, axis=0, mode='raise').copy('F')
        n_upper = self.level_populations.values.take(self.atom_data.lines_upper2level_idx, axis=0, mode='raise').copy('F')
        meta_stable_upper = self.atom_data.levels.metastable.values.take(self.atom_data.lines_upper2level_idx, axis=0, mode='raise')
        g_lower = self.atom_data.levels.g.values.take(self.atom_data.lines_lower2level_idx, axis=0, mode='raise')
        g_upper = self.atom_data.levels.g.values.take(self.atom_data.lines_upper2level_idx, axis=0, mode='raise')


        self.stimulated_emission_factor = 1 - ((g_lower[np.newaxis].T * n_upper) / (g_upper[np.newaxis].T * n_lower))
        # getting rid of the obvious culprits
        self.stimulated_emission_factor[n_lower == 0.0] = 0.0
        self.stimulated_emission_factor[np.isneginf(self.stimulated_emission_factor)] = 0.0
        self.stimulated_emission_factor[meta_stable_upper[np.newaxis].T & (self.stimulated_emission_factor < 0)] = 0.0

        if self.nlte_config is not None and self.nlte_config.species != []:
            nlte_lines_mask = np.zeros(self.stimulated_emission_factor.shape[0]).astype(bool)
            for species in self.nlte_config.species:
                nlte_lines_mask |= (self.atom_data.lines.atomic_number == species[0]) & \
                                   (self.atom_data.lines.ion_number == species[1])
            self.stimulated_emission_factor[(self.stimulated_emission_factor < 0) & nlte_lines_mask[np.newaxis].T] = 0.0


        tau_sobolevs = sobolev_coefficient * f_lu[np.newaxis].T * wavelength[np.newaxis].T * self.time_explosion * \
                       n_lower * self.stimulated_emission_factor
        return pd.DataFrame(tau_sobolevs, index=self.atom_data.lines.index, columns=np.arange(len(self.t_rads)))




    def calculate_transition_probabilities(self):
        """
            Updating the Macro Atom computations
        """

        macro_atom_data = self.atom_data.macro_atom_data
        if not hasattr(self, 'beta_sobolevs'):
            self.beta_sobolevs = np.zeros_like(self.tau_sobolevs.values)

        if not self.beta_sobolevs_precalculated:
            macro_atom.calculate_beta_sobolev(self.tau_sobolevs.values.ravel(order='F'),
                                          self.beta_sobolevs.ravel(order='F'))

        transition_probabilities = (macro_atom_data.transition_probability.values[np.newaxis].T *
                                    self.beta_sobolevs.take(self.atom_data.macro_atom_data.lines_idx.values.astype(int),
                                                            axis=0, mode='raise')).copy('F')
        transition_up_filter = (macro_atom_data.transition_type == 1).values
        macro_atom_transition_up_filter = macro_atom_data.lines_idx.values[transition_up_filter]
        j_blues = self.j_blues.values.take(macro_atom_transition_up_filter, axis=0, mode='raise')
        macro_stimulated_emission = self.stimulated_emission_factor.take(macro_atom_transition_up_filter, axis=0, mode='raise')
        transition_probabilities[transition_up_filter] *= j_blues * macro_stimulated_emission
        #Normalizing the probabilities
        block_references = np.hstack((self.atom_data.macro_atom_references.block_references,
                                      len(macro_atom_data)))
        macro_atom.normalize_transition_probabilities(transition_probabilities, block_references)
        return pd.DataFrame(transition_probabilities, index=macro_atom_data.transition_line_id,
                     columns=self.tau_sobolevs.columns)
