import pandas as pd
import numpy as np

from tardis.continuum.base import ContinuumProcess
from tardis import macro_atom


class BaseTransitionProbabilities(ContinuumProcess):
    def __init__(self, input_data, **kwargs):
        super(BaseTransitionProbabilities, self).__init__(input_data)
        if not kwargs:
            raise ValueError
        self.dataframe = self._prepare_transition_probabilities(**kwargs)
        self.data_array = None
        self.block_references = None
        self._set_montecarlo_data()

        def _prepare_transition_probabilities(self, **kwargs):
            pass

        def _set_montecarlo_data(self):
            pass


class TransitionProbabilities(BaseTransitionProbabilities):
    def __init__(self, input_data, add_probabilities=True, **kwargs):
        self.add_probabilities = add_probabilities
        self.destination_level_id = None
        self.transition_type = None
        self.transition_line_id = None
        super(TransitionProbabilities, self).__init__(input_data, **kwargs)

    def _prepare_transition_probabilities(self, **kwargs):
        mode2append = {'up': 'internal_up', 'down': 'internal_down', 'continuum': 'continuum'}
        transition_probabilities_dict = {'internal_up': [], 'internal_down': [], 'radiative_deactivation': [],
                                         'collisional_deactivation': [], 'continuum': []}
        for name, process in kwargs.iteritems():
            internal_jump = process.internal_jump_probabilities
            if process.macro_atom_transitions == 'down':
                deactivation = process.deactivation_probabilities
                if process.name == 'collisional_deexcitation':
                    transition_probabilities_dict['collisional_deactivation'].append(deactivation)
                elif process.name == 'radiative_deexcitation':
                    transition_probabilities_dict['radiative_deactivation'].append(deactivation)
                else:
                    raise ValueError

            transition_probabilities_dict[mode2append[process.macro_atom_transitions]].append(internal_jump)

        if self.add_probabilities == True:
            transition_probabilities_dict = self._add_probabilities(transition_probabilities_dict)

        transition_probabilities_dict = self._set_references(transition_probabilities_dict)

        transition_probabilities = self._combine_transition_probabilities(transition_probabilities_dict)
        return transition_probabilities


    def _set_references(self, transition_probabilities_dict):
        transition_type2idx = {'internal_up': 1, 'internal_down': 0, 'continuum': 2, 'radiative_deactivation': -1,
                               'collisional_deactivation': -2}
        for transition_type, probabilities_list in transition_probabilities_dict.iteritems():
            for probabilities in probabilities_list:
                ones = self.ones(probabilities)
                probabilities.insert(0, 'transition_type', transition_type2idx[transition_type] * ones)
                if not 'lines_idx' in probabilities.columns:
                    probabilities.insert(1, 'lines_idx', -1 * ones)
                if transition_type == 'continuum':
                    import ipdb;

                    ipdb.set_trace()
                    multi_index = self._get_ion_prob_index(probabilities.index)
                    probabilities.set_index(multi_index, inplace=True)
        return transition_probabilities_dict


    def _combine_transition_probabilities(self, transition_probabilities_dict):
        transition_probabilities_list = []
        for transition_type, probabilities in transition_probabilities_dict.iteritems():
            transition_probabilities_list.extend(probabilities)
        combined_probabilities = pd.concat(transition_probabilities_list)
        combined_probabilities.sortlevel(sort_remaining=False, inplace=True)
        self.block_references = self._get_block_references(combined_probabilities)
        combined_probabilities = self._normalize_transition_probabilities(
            combined_probabilities, no_ref_columns=2)
        return combined_probabilities

    def _set_montecarlo_data(self):
        self.block_references = self._get_block_references(self.dataframe)
        self.data_array = np.ascontiguousarray(self.dataframe.ix[:, 2:].values)
        self.transition_type = self.dataframe['transition_type'].values.astype(np.int64)
        self.transition_line_id = self.dataframe['lines_idx'].values.astype(np.int64)
        self.destination_level_id = self.dataframe.index.get_level_values(1).values

    def _add_probabilities(self, transition_probabilities_dict):
        for key in ['internal_up', 'internal_down', 'continuum']:
            probabilities2add = transition_probabilities_dict[key]
            if len(probabilities2add) == 2:
                added_probabilities = probabilities2add[0].add(probabilities2add[1], fill_value=0.)
            elif probabilities2add == 1:
                added_probabilities = probabilities2add[0]
            else:
                raise ValueError

            transition_probabilities_dict[key] = [added_probabilities]

        return transition_probabilities_dict


class RecombinationTransitionProbabilities(BaseTransitionProbabilities):
    def __init__(self, input_data, **kwargs):
        super(RecombinationTransitionProbabilities, self).__init__(input_data, **kwargs)

    def _prepare_transition_probabilities(self, **kwargs):
        transition_probabilities = []
        transition_probabilities.append(self._calculate_internal_jump_probabilities(**kwargs))
        transition_probabilities.extend(self._calculate_deactivation_probabilities(**kwargs))

        transition_probabilities = pd.concat(transition_probabilities)

        transition_probabilities = self._normalize_transition_probabilities(
            transition_probabilities, no_ref_columns=3)
        return transition_probabilities

    def _set_montecarlo_data(self):
        self.data_array = self._get_contiguous_array(self.dataframe.ix[:, 3:])
        self.data_array_nd = self.data_array.shape[1]
        self.block_references = self._get_block_references(self.dataframe)

    def _calculate_internal_jump_probabilities(self, **kwargs):
        if len(kwargs) == 2:
            trans_prob = \
                kwargs['radiative_recombination'].internal_jump_probabilities.add(
                    kwargs['collisional_recombination'].internal_jump_probabilities)
        else:
            trans_prob = kwargs.values()[0].internal_jump_probabilities

        destination_level_idx = self._get_level_idx(trans_prob.index)

        trans_prob.insert(0, 'destination_level_idx', destination_level_idx)
        trans_prob.insert(1, 'continuum_edge_idx', - 1 * self.ones(trans_prob))
        trans_prob.insert(2, 'transition_type', 0 * self.ones(trans_prob))
        return trans_prob

    def _calculate_deactivation_probabilities(self, **kwargs):
        deactivation_probabilities = []
        for name, process in kwargs.iteritems():
            trans_prob = process.deactivation_probabilities
            self._set_deactivation_references(trans_prob, process.name)
            deactivation_probabilities.append(trans_prob)
        return deactivation_probabilities

    def _set_deactivation_references(self, deact_prob, process_name):
        destination_level_idx = - 1 * self.ones(deact_prob)
        if process_name == 'collisional_recombination':
            continuum_edge_idx = - 1 * self.ones(deact_prob)
            transition_type = -2 * self.ones(deact_prob)
        else:
            continuum_edge_idx = self._get_continuum_edge_idx(deact_prob.index)
            transition_type = -3 * self.ones(deact_prob)

        deact_prob.insert(0, 'destination_level_idx', destination_level_idx)
        deact_prob.insert(1, 'continuum_edge_idx', continuum_edge_idx)
        deact_prob.insert(2, 'transition_type', transition_type)
        return deact_prob


class LegacyTransitionProbabilities(ContinuumProcess):
    def __init__(self, input_data, radiative_prob, coll_int_jump_up_prob, coll_int_jump_down_prob,
                 coll_deact_prob, coll_ion_prob, rad_ion_prob):
        super(LegacyTransitionProbabilities, self).__init__(input_data)

        self.destination_level_id = None
        self.transition_type = None
        self.transition_line_id = None

        self.radiative_probabilities = \
            self._prepare_radiative_probabilities(radiative_prob)
        self.collisional_deactivation_probabilities = \
            self._prepare_collisional_deactivation_probabilities(coll_deact_prob)
        self.ionization_probabilities = self._prepare_ionization_probabilities(coll_ion_prob, rad_ion_prob)
        # This is mostly useful when using Van Regemorter
        self.added_probabilities = self._add_internal_jump_probabilities(coll_int_jump_up_prob=coll_int_jump_up_prob,
                                                                         coll_int_jump_down_prob=coll_int_jump_down_prob)
        self.dataframe = self._combine_transition_probabilities(
            self.added_probabilities, self.collisional_deactivation_probabilities, self.ionization_probabilities)

        self._set_montecarlo_data()


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

    def _prepare_ionization_probabilities(self, coll_ion_prob, rad_ion_prob):
        # WARNING: destination level id is continuum id; the value itself is not unique
        ion_prob_prep = coll_ion_prob.add(rad_ion_prob)
        ion_prob_prep.insert(0, 'transition_type', 2 * np.ones(ion_prob_prep.values.shape[0]))
        ion_prob_prep.insert(1, 'lines_idx', -1 * np.ones(ion_prob_prep.values.shape[0]))
        multi_index = self._get_ion_prob_index(ion_prob_prep.index)
        ion_prob_prep.set_index(multi_index, inplace=True)
        return ion_prob_prep

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
        self.block_references = self._get_block_references(combined_probabilities)
        macro_atom.normalize_transition_probabilities(combined_probabilities.ix[:, 2:].values, self.block_references)
        return combined_probabilities

    def _set_montecarlo_data(self):
        self.block_references = self._get_block_references(self.dataframe)
        self.data_array = np.ascontiguousarray(self.dataframe.ix[:, 2:].values)
        self.transition_type = self.dataframe['transition_type'].values.astype(np.int64)
        self.transition_line_id = self.dataframe['lines_idx'].values.astype(np.int64)
        self.destination_level_id = self.dataframe.index.get_level_values(1).values

    def _get_source_level_idx(self):
        macro_atom_data = self.macro_atom_data
        source_level_index = pd.MultiIndex.from_arrays([macro_atom_data['atomic_number'], macro_atom_data['ion_number'],
                                                        macro_atom_data['source_level_number']])
        return self._get_level_idx(source_level_index)

    def _get_ion_prob_index(self, level_lower_index):
        source_level_idx = self._get_level_idx(level_lower_index)
        destination_level_idx = self._get_continuum_idx(level_lower_index)
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        return tmp_multi_index
