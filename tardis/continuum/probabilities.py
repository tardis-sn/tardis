from tardis.continuum.base import ContinuumProcess
import pandas as pd
import numpy as np
from tardis import macro_atom


class TransitionProbabilitiesCombined(ContinuumProcess):
    def __init__(self, input_data, radiative_prob, coll_int_jump_up_prob, coll_int_jump_down_prob,
                 coll_deact_prob, ion_prob):
        super(TransitionProbabilitiesCombined, self).__init__(input_data)

        self.combined_transition_probabilities = None
        self.block_references = None
        self.destination_level_id = None
        self.transition_probabilities_array = None
        self.radiative_probabilities = \
            self._prepare_radiative_probabilities(radiative_prob)
        self.collisional_deactivation_probabilities = \
            self._prepare_collisional_deactivation_probabilities(coll_deact_prob)
        self.ionization_probabilities = self._prepare_ionization_probabilities(ion_prob)
        # This is mostly useful when using Van Regemorter
        self.added_probabilities = self._add_internal_jump_probabilities(coll_int_jump_up_prob=coll_int_jump_up_prob,
                                                                         coll_int_jump_down_prob=coll_int_jump_down_prob)
        self._combine_transition_probabilities(self.added_probabilities, self.collisional_deactivation_probabilities,
                                               self.ionization_probabilities)

    def _prepare_radiative_probabilities(self, radiative_prob):
        source_level_idx = self._get_source_level_idx()
        destination_level_idx = self.macro_atom_data.destination_level_idx.values
        new_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                              names=['source_level_idx', 'destination_level_idx'])
        radiative_prob_prep = radiative_prob.set_index(new_index)
        radiative_prob_prep.insert(0, 'transition_type', self.macro_atom_data.transition_type.values)
        radiative_prob_prep.insert(1, 'lines_idx', self.macro_atom_data['lines_idx'].values)
        return radiative_prob_prep

    def _prepare_collisional_deactivation_probabilities(self, coll_deact_prob):
        coll_deact_prob_prep = coll_deact_prob.copy()
        coll_deact_prob_prep.insert(0, 'transition_type', -2 * np.ones(coll_deact_prob.values.shape[0]))
        coll_deact_prob_prep.insert(1, 'lines_idx', -1 * np.ones(coll_deact_prob.values.shape[0]))
        return coll_deact_prob_prep

    def _prepare_ionization_probabilities(self, ion_prob):
        # WARNING: destination level id is continuum id; the value itself is not unique
        ion_prob_prep = ion_prob.copy()
        ion_prob_prep.insert(0, 'transition_type', 2 * np.ones(ion_prob.values.shape[0]))
        ion_prob_prep.insert(1, 'lines_idx', -1 * np.ones(ion_prob.values.shape[0]))
        multi_index = self._get_ion_prob_index(ion_prob_prep.index)
        ion_prob_prep.set_index(multi_index, inplace=True)
        return ion_prob_prep

    def _get_source_level_idx(self):
        import ipdb;

        ipdb.set_trace()
        macro_atom_data = self.macro_atom_data
        source_level_index = pd.MultiIndex.from_arrays([macro_atom_data['atomic_number'], macro_atom_data['ion_number'],
                                                        macro_atom_data['source_level_number']])
        return self._get_level_idx(source_level_index)

    def _add_internal_jump_probabilities(self, coll_int_jump_up_prob, coll_int_jump_down_prob):
        added_probabilities = self.radiative_probabilities.copy()
        transition_up_filter = self.transition_up_filter
        transition_down_filter = self.transition_down_filter

        added_probabilities[transition_up_filter] = \
            added_probabilities[transition_up_filter].combineAdd(coll_int_jump_up_prob)
        added_probabilities[transition_down_filter] = \
            added_probabilities[transition_down_filter].combineAdd(coll_int_jump_down_prob)
        return added_probabilities

    def _combine_transition_probabilities(self, *args):
        combined_probabilities = pd.concat(args)
        combined_probabilities.sortlevel(sort_remaining=False, inplace=True)
        self.block_references = self._get_new_block_references(combined_probabilities)
        macro_atom.normalize_transition_probabilities(combined_probabilities.ix[:, 2:].values, self.block_references)
        self.combined_transition_probabilities = combined_probabilities
        self.destination_level_id = combined_probabilities.index.get_level_values(1).values
        # ? transpose
        self.transition_probabilities_array = np.ascontiguousarray(
            self.combined_transition_probabilities.ix[:, 2:].values)
        self.transition_type = self.combined_transition_probabilities['transition_type'].values.astype(np.int64)
        self.transition_line_id = self.combined_transition_probabilities['lines_idx'].values.astype(np.int64)

    def _get_new_block_references(self, combined_probabilities):
        block_references = combined_probabilities[0].groupby(level=0).count().cumsum().values
        block_references = np.hstack([[0], block_references])
        return block_references

    def _get_ion_prob_index(self, level_lower_index):
        source_level_idx = self._get_level_idx(level_lower_index)
        destination_level_idx = self._get_continuum_idx(level_lower_index)
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        return tmp_multi_index