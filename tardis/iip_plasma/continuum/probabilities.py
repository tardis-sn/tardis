import numpy as np
import pandas as pd

from tardis.iip_plasma.continuum.base import ContinuumProcess


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
        mode2append = {
            "up": "internal_up",
            "down": "internal_down",
            "continuum": "continuum",
        }
        transition_probabilities_dict = {
            "internal_up": [],
            "internal_down": [],
            "radiative_deactivation": [],
            "collisional_deactivation": [],
            "continuum": [],
        }
        for name, process in kwargs.iteritems():
            internal_jump = process.internal_jump_probabilities
            if process.macro_atom_transitions == "down":
                deactivation = process.deactivation_probabilities
                if process.name == "collisional_deexcitation":
                    transition_probabilities_dict[
                        "collisional_deactivation"
                    ].append(deactivation)
                elif process.name == "radiative_deexcitation":
                    transition_probabilities_dict[
                        "radiative_deactivation"
                    ].append(deactivation)
                else:
                    raise ValueError

            transition_probabilities_dict[
                mode2append[process.macro_atom_transitions]
            ].append(internal_jump)

        if self.add_probabilities == True:
            transition_probabilities_dict = self._add_probabilities(
                transition_probabilities_dict
            )

        transition_probabilities_dict = self._set_references(
            transition_probabilities_dict
        )

        transition_probabilities = self._combine_transition_probabilities(
            transition_probabilities_dict
        )
        return transition_probabilities

    def _set_references(self, transition_probabilities_dict):
        transition_type2idx = {
            "internal_up": 1,
            "internal_down": 0,
            "continuum": 2,
            "radiative_deactivation": -1,
            "collisional_deactivation": -2,
        }
        for (
            transition_type,
            probabilities_list,
        ) in transition_probabilities_dict.iteritems():
            for probabilities in probabilities_list:
                ones = self.ones(probabilities)
                probabilities.insert(
                    0,
                    "transition_type",
                    transition_type2idx[transition_type] * ones,
                )
                if "lines_idx" not in probabilities.columns:
                    probabilities.insert(1, "lines_idx", -1 * ones)
                if transition_type == "continuum":
                    multi_index = self._get_ion_prob_index(probabilities.index)
                    probabilities.set_index(multi_index, inplace=True)
        return transition_probabilities_dict

    def _combine_transition_probabilities(self, transition_probabilities_dict):
        transition_probabilities_list = []
        for (
            transition_type,
            probabilities,
        ) in transition_probabilities_dict.iteritems():
            transition_probabilities_list.extend(probabilities)
        combined_probabilities = pd.concat(transition_probabilities_list)
        combined_probabilities.sortlevel(sort_remaining=False, inplace=True)
        self.block_references = self._get_block_references(
            combined_probabilities
        )
        combined_probabilities = self._normalize_transition_probabilities(
            combined_probabilities, no_ref_columns=2
        )
        return combined_probabilities

    def _set_montecarlo_data(self):
        self.block_references = self._get_block_references(self.dataframe)
        self.data_array = np.ascontiguousarray(
            self.dataframe.iloc[:, 2:].values
        )
        self.transition_type = self.dataframe["transition_type"].values.astype(
            np.int64
        )
        self.transition_line_id = self.dataframe["lines_idx"].values.astype(
            np.int64
        )
        self.destination_level_id = self.dataframe.index.get_level_values(
            1
        ).values

    def _add_probabilities(self, transition_probabilities_dict):
        for key in ["internal_up", "internal_down", "continuum"]:
            probabilities2add = transition_probabilities_dict[key]
            if len(probabilities2add) == 2:
                added_probabilities = probabilities2add[0].add(
                    probabilities2add[1], fill_value=0.0
                )
            elif probabilities2add == 1:
                added_probabilities = probabilities2add[0]
            else:
                raise ValueError

            transition_probabilities_dict[key] = [added_probabilities]

        return transition_probabilities_dict

    def _get_block_references(self, probabilities):
        block_references_series = (
            probabilities[0].groupby(level=0).count().cumsum()
        )
        # TODO: Check if this works after the bug fix
        no_of_levels = self.input.atom_data.levels.shape[0]
        block_references = np.zeros(no_of_levels + 1, dtype=np.int64)
        # Necessary because metastable levels do not have transition probabilities
        block_references[block_references_series.index + 1] = (
            block_references_series.values
        )
        last_nonzero_entry = 0
        for i in range(1, len(block_references)):
            if block_references[i] == 0:
                block_references[i] = last_nonzero_entry
            else:
                last_nonzero_entry = block_references[i]
        return block_references


class RecombinationTransitionProbabilities(BaseTransitionProbabilities):
    def __init__(self, input_data, **kwargs):
        super(RecombinationTransitionProbabilities, self).__init__(
            input_data, **kwargs
        )

    def _prepare_transition_probabilities(self, **kwargs):
        transition_probabilities = []
        transition_probabilities.append(
            self._calculate_internal_jump_probabilities(**kwargs)
        )
        transition_probabilities.extend(
            self._calculate_deactivation_probabilities(**kwargs)
        )

        transition_probabilities = pd.concat(transition_probabilities)

        transition_probabilities = self._normalize_transition_probabilities(
            transition_probabilities, no_ref_columns=3
        )
        return transition_probabilities

    def _set_montecarlo_data(self):
        self.data_array = self._get_contiguous_array(self.dataframe.iloc[:, 3:])
        self.data_array_nd = self.data_array.shape[1]
        self.block_references = self._get_block_references(self.dataframe)

    def _calculate_internal_jump_probabilities(self, **kwargs):
        if len(kwargs) == 2:
            trans_prob = kwargs[
                "radiative_recombination"
            ].internal_jump_probabilities.add(
                kwargs["collisional_recombination"].internal_jump_probabilities
            )
        else:
            trans_prob = kwargs.values()[0].internal_jump_probabilities

        destination_level_idx = self._get_level_idx(trans_prob.index)

        trans_prob.insert(0, "destination_level_idx", destination_level_idx)
        trans_prob.insert(1, "continuum_edge_idx", -1 * self.ones(trans_prob))
        trans_prob.insert(2, "transition_type", 0 * self.ones(trans_prob))
        return trans_prob

    def _calculate_deactivation_probabilities(self, **kwargs):
        deactivation_probabilities = []
        for name, process in kwargs.iteritems():
            trans_prob = process.deactivation_probabilities
            self._set_deactivation_references(trans_prob, process.name)
            deactivation_probabilities.append(trans_prob)
        return deactivation_probabilities

    def _set_deactivation_references(self, deact_prob, process_name):
        destination_level_idx = -1 * self.ones(deact_prob)
        if process_name == "collisional_recombination":
            continuum_edge_idx = -1 * self.ones(deact_prob)
            transition_type = -2 * self.ones(deact_prob)
        else:
            continuum_edge_idx = self._get_continuum_edge_idx(deact_prob.index)
            transition_type = -3 * self.ones(deact_prob)

        deact_prob.insert(0, "destination_level_idx", destination_level_idx)
        deact_prob.insert(1, "continuum_edge_idx", continuum_edge_idx)
        deact_prob.insert(2, "transition_type", transition_type)
        return deact_prob

    def _get_block_references(self, probabilities):
        block_references = (
            probabilities[0].groupby(level=0).count().cumsum().values
        )
        block_references = np.hstack([[0], block_references])
        return block_references
