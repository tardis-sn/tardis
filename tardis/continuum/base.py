import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as units


def data_type_selection(data_getter):
    def wrapped_data_getter(*args, **kwargs):
        output = data_getter(*args)
        if not kwargs or kwargs['dtype'] == 'array':
            return output.values
        elif kwargs['dtype'] == 'dataframe':
            return output
        else:
            raise AttributeError

    return wrapped_data_getter


class ContinuumProcess(object):
    def __init__(self, input_data):
        self.input = input_data

    @data_type_selection
    def _get_level_energy(self, multi_index):
        return self.input.levels.loc[multi_index, 'energy']

    @data_type_selection
    def _get_lte_level_pop(self, multi_index):
        return self.input.lte_level_pop.loc[multi_index]

    @data_type_selection
    def _get_level_pop(self, multi_index):
        return self.input.level_pop.loc[multi_index]

    def _get_level_idx(self, multi_index):
        return self.input.macro_atom_references.loc[multi_index, 'references_idx'].values

    def _get_continuum_idx(self, multi_index_full):
        ion_number_index = self._get_ion_multi_index(multi_index_full, next_higher=False)
        return self.input.continuum_references.loc[ion_number_index, 'references_idx'].values

    @staticmethod
    def _get_ion_multi_index(multi_index_full, next_higher=True):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1)
        if next_higher is True:
            ion_number += 1
        return pd.MultiIndex.from_arrays([atomic_number, ion_number])

    @data_type_selection
    def _get_ion_number_density(self, multi_index_full):
        ion_number_index = self._get_ion_multi_index(multi_index_full)
        return self.ion_number_density.loc[ion_number_index]

    @data_type_selection
    def _get_lte_ion_number_density(self, multi_index_full):
        ion_number_index = self._get_ion_multi_index(multi_index_full)
        # WARNING: In the future we have to return the lte ion number density
        return self.ion_number_density.loc[ion_number_index]

    @data_type_selection
    def _get_ionization_energy(self, multi_index_full):
        ion_number_index = self._get_ion_multi_index(multi_index_full)
        return self.input.ionization_energies.loc[ion_number_index, 'ionization_energy']

    def _get_continuum_edge_idx(self, multi_index):
        return self.input.continuum_data.set_index(['atomic_number', 'ion_number',
                                                    'level_number_lower']).loc[multi_index, 'continuum_edge_idx']

    def _get_block_references(self, probabilities):
        block_references = probabilities[0].groupby(level=0).count().cumsum().values
        block_references = np.hstack([[0], block_references])
        return block_references

    @staticmethod
    def _normalize_transition_probabilities(dataframe, no_ref_columns=0):
        normalization_fct = lambda x: (x / x.sum())
        normalized_dataframe = dataframe.ix[:, no_ref_columns:].groupby(level=0).transform(normalization_fct)
        normalized_dataframe = pd.concat([dataframe.ix[:, :no_ref_columns], normalized_dataframe], axis=1,
                                         join_axes=[normalized_dataframe.index])
        normalized_dataframe = normalized_dataframe.fillna(0.0)
        return normalized_dataframe

    def _get_level_multi_index(self, level_idx):
        tmp = self.input.macro_atom_references_by_idx.iloc[level_idx]
        level_multi_index = pd.MultiIndex.from_arrays([tmp['atomic_number'].values, tmp['ion_number'].values,
                                                       tmp['source_level_number'].values])
        return level_multi_index

    def _set_ionization_rates_index(self, probabilities):
        # WARNING: destination level id is continuum id; the value itself is not unique
        multi_index = self._get_ion_prob_index(probabilities.index)
        probabilities.set_index(multi_index, inplace=True)

    def _get_ion_prob_index(self, level_lower_index):
        source_level_idx = self._get_level_idx(level_lower_index)
        destination_level_idx = self._get_continuum_idx(level_lower_index)
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        return tmp_multi_index

    @property
    def electron_densities(self):
        return self.input.electron_densities

    @property
    def t_rads(self):
        return self.input.t_rads

    @property
    def ws(self):
        return self.input.ws

    @property
    def t_electrons(self):
        return self.input.t_electrons

    @property
    def ion_number_density(self):
        return self.input.ion_number_density

    @property
    def transition_up_filter(self):
        return (self.input.macro_atom_data.transition_type == 1).values

    @property
    def transition_down_filter(self):
        return (self.input.macro_atom_data.transition_type == 0).values

    @property
    def transition_deactivation_filter(self):
        return (self.input.macro_atom_data.transition_type == -1).values

    @property
    def macro_atom_data(self):
        return self.input.macro_atom_data

    @property
    def nu_i(self):
        return self.input.nu_i

    @property
    def photoionization_data(self):
        return self.input.photoionization_data

    @property
    def no_of_shells(self):
        return len(self.input.t_rads)

    # Helper functions for physical calculations
    def _calculate_u0s(self, nu):
        u0s = nu[np.newaxis].T / self.t_electrons * (const.h.cgs.value / const.k_B.cgs.value)
        return u0s

    # Data preparation
    def _get_contiguous_array(self, dataframe):
        return np.ascontiguousarray(dataframe.values.transpose())

    @staticmethod
    def ones(dataframe, dtype=np.int64):
        return np.ones(dataframe.shape[0], dtype)


class TransitionProbabilitiesMixin(object):
    @property
    def internal_jump_probabilities(self):
        level_lower_energy = self.level_lower_energy * units.erg.to(units.eV)
        return self.rate_coefficient.multiply(level_lower_energy, axis=0)  # / cconst.c_einstein

    @property
    def deactivation_probabilities(self):
        energy_difference = (self.level_upper_energy - self.level_lower_energy) * units.erg.to(units.eV)
        return self.rate_coefficient.multiply(energy_difference, axis=0)  # / cconst.c_einstein


class PhysicalContinuumProcess(ContinuumProcess, TransitionProbabilitiesMixin):
    name = None
    cooling = True
    macro_atom_transitions = None

    def __init__(self, input_data, **kwargs):
        super(PhysicalContinuumProcess, self).__init__(input_data)
        self.rate_coefficient = self._calculate_rate_coefficient(**kwargs)
        if self.cooling is True:
            self.cooling_rate = self._calculate_cooling_rate(**kwargs)

    def _calculate_rate_coefficient(self, **kwargs):
        return None

    def _calculate_cooling_rate(self, **kwargs):
        raise NotImplementedError


class BoundFreeEnergyMixIn(object):
    @property
    def level_lower_energy(self):
        return self._get_level_energy(self.rate_coefficient.index)

    @property
    def level_upper_energy(self):
        return self._get_ionization_energy(self.rate_coefficient.index)


class InverseProcess(ContinuumProcess, TransitionProbabilitiesMixin):
    name = None
    name_of_inverse_process = None
    cooling = False
    macro_atom_transitions = None

    def __init__(self, input_data, rate_coefficient, inverse_process):
        self.input = input_data
        self.rate_coefficient = rate_coefficient
        self.inverse_process = inverse_process

    @classmethod
    def from_inverse_process(cls, inverse_process):
        rate_coefficient = cls._calculate_inverse_rate(inverse_process)
        return cls(inverse_process.input, rate_coefficient, inverse_process)

    @classmethod
    def _calculate_inverse_rate(cls, inverse_process):
        raise NotImplementedError

    @property
    def level_lower_energy(self):
        return self.inverse_process.level_lower_energy

    @property
    def level_upper_energy(self):
        return self.inverse_process.level_upper_energy

        # @property
        # def transition_probabilities_int_down(self):
        #    return self.rate_coefficient.multiply(self.inverse_process.energy_lower, axis=0) / self.c_einstein

        #@property
        #def transition_probabilities_deactivation(self):
        #    return self.rate_coefficient.multiply(self.inverse_process.energy_difference, axis=0) / self.c_einstein
