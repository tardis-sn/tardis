import numpy as np
import logging
from pandas import DataFrame
from astropy import units
from astropy import constants as const
import pandas as pd
from scipy.integrate import simps
from tardis import macro_atom
from tardis.util import intensity_black_body
from scipy.special import expn

logger = logging.getLogger(__name__)


class IonizationRates(object):
    def __init__(self, photoionization_data, ws, t_rads, electron_densities, level_pop, continuum_data):
        self.photoionization_data = photoionization_data
        self.nu_i = self._get_nu_i()
        self.continuum_data = continuum_data
        self.ws = ws
        self.t_rads = t_rads
        self.level_pop = level_pop
        self.electron_densities = electron_densities
        self.corrected_photoionization_coefficient = self._calculate_corrected_photoion_coeff()
        self.collisional_ionization_coefficient = self._calculate_collion_coeff()
        self.coll_ionization_cooling_rate = self._calculate_coll_ionization_cooling_rate()

    def _calculate_corrected_photoion_coeff(self):
        j_nus = self._calculate_j_nus()
        stimulated_emission_correction = self._calculate_stimulated_emission_correction()
        corrected_photoion_coeff = j_nus.multiply(4. * np.pi * self.photoionization_data['x_sect'] /
                                                  self.photoionization_data['nu'] / const.h.cgs.value, axis=0)
        corrected_photoion_coeff = corrected_photoion_coeff.multiply(stimulated_emission_correction)
        corrected_photoion_coeff.insert(0, 'nu', self.photoionization_data['nu'])
        corrected_photoion_coeff = corrected_photoion_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(self.t_rads)):
            tmp[i] = corrected_photoion_coeff.apply(lambda sub: simps(sub[i], sub['nu'], even='first'))
        corrected_photoion_coeff = pd.DataFrame(tmp)
        return corrected_photoion_coeff

    def _calculate_j_nus(self):
        nus = self.photoionization_data['nu'].values
        j_nus = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads.value)
        return pd.DataFrame(j_nus, index=self.photoionization_data.index, columns=np.arange(len(self.t_rads)))

    def _calculate_collion_coeff(self):
        collion_coeff = self.photoionization_data.groupby(level=[0, 1, 2]).first()
        factor = self._calculate_factor(collion_coeff['nu'].values, index=collion_coeff.index)
        collion_coeff = 1.55e13 * collion_coeff['x_sect']
        ion_number = collion_coeff.index.get_level_values(1).values
        collion_coeff.ix[ion_number == 0] *= 0.1
        collion_coeff.ix[ion_number == 1] *= 0.2
        collion_coeff.ix[ion_number >= 2] *= 0.3
        collion_coeff = factor.multiply(collion_coeff, axis=0)
        collion_coeff = collion_coeff.multiply(self.electron_densities / np.sqrt(0.9 * self.t_rads.value), axis=1)
        return collion_coeff

    def _calculate_boltzmann_factor(self, nu):
        u0s = self._calculate_u0s(nu)
        return np.exp(-u0s)

    def _calculate_u0s(self, nu):
        # WARNING: This might be unsafe for future implementations.
        t_electrons = 0.9 * self.t_rads.value
        u0s = nu[np.newaxis].T / t_electrons * (const.h.cgs.value / const.k_B.cgs.value)
        return u0s

    def _calculate_factor(self, nu_ijk, index):
        u0s = self._calculate_u0s(nu_ijk)
        factor = 1 / u0s * np.exp(-u0s)
        factor = pd.DataFrame(factor, index=index, columns=np.arange(len(self.t_rads)))
        return factor

    def _calculate_stimulated_emission_correction(self):
        nu = self.photoionization_data['nu'].values
        boltzmann_factor = self._calculate_boltzmann_factor(nu)
        # TODO: generalize
        lte_nonlte_level_pop_ratio = 1.
        correction_factor = (1. - lte_nonlte_level_pop_ratio * boltzmann_factor)
        return correction_factor

    def _calculate_coll_ionization_cooling_rate(self):
        level_pop = self.level_pop.loc[self.collisional_ionization_coefficient.index]
        coll_ion_cooling_rate = self.collisional_ionization_coefficient.multiply(level_pop)
        coll_ion_cooling_rate = coll_ion_cooling_rate.multiply(self.nu_i * const.h.cgs.value, axis=0)
        continuum_idx = self._get_continuum_idx(coll_ion_cooling_rate.index)
        coll_ion_cooling_rate.set_index(continuum_idx, inplace=True)
        return coll_ion_cooling_rate

    def _get_continuum_idx(self, multi_index_full):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1)
        ion_number_index = pd.MultiIndex.from_arrays([atomic_number, ion_number])
        return self.continuum_data.continuum_references.loc[ion_number_index, 'references_idx']

    def _get_nu_i(self):
        return self.photoionization_data.groupby(level=[0, 1, 2]).first().nu.values

    @property
    def transition_probabilities2cont(self):
        c_einstein = (4. * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)).value
        trans_prob2cont = \
            self.corrected_photoionization_coefficient.add(self.collisional_ionization_coefficient) / c_einstein
        # TODO check that this is safe under all circumstances
        trans_prob2cont = trans_prob2cont.multiply(self.nu_i, axis=0) * const.h.cgs.value * units.erg.to(units.eV)
        return trans_prob2cont


class CollisionalRates(object):
    def __init__(self, lines, t_electrons, electron_densities, lte_level_pop,
                 levels, macro_atom_references, mode='Van Regemorter'):
        self.t_electrons = t_electrons
        self.n_e = electron_densities
        self.mode = mode
        if mode == 'Van Regemorter':
            self.level_lower_index = self._get_level_lower_index(lines)
            self.level_upper_index = self._get_level_upper_index(lines)
            self.coll_excitation_coeff = self.calculate_coll_excitation_coeff_regemorter(lines=lines,
                                                                                         macro_atom_references=macro_atom_references)
        else:
            raise NotImplementedError

        # TODO: replace lte level pop
        self.excitation_cooling_rate = self._calculate_excitation_cooling_rate(level_pop=lte_level_pop, levels=levels,
                                                           coll_excitation_coeff=self.coll_excitation_coeff)

        self.coll_deexcitation_coeff = self.calculate_coll_deexcitation_coeff(
            coll_excitation_coeff=self.coll_excitation_coeff, lte_level_pop=lte_level_pop)
        self.transition_prob_int_up, self.transition_prob_int_down, self.transition_prob_deact = \
            self.transition_probabilities(coll_excitation_coeff=self.coll_excitation_coeff,
                                          coll_deexcitation_coeff=self.coll_deexcitation_coeff, levels=levels)

    def calculate_coll_excitation_coeff_regemorter(self, lines, macro_atom_references):
        c_0 = 5.465e-11
        I_H = 13.598433770784 * units.eV.to(units.erg)
        coll_excitation_coeff = pd.DataFrame(14.5 * c_0 * lines.f_lu.values *
                                             (I_H / (const.h.cgs.value * lines.nu.values)) ** 2)
        coll_excitation_coeff = coll_excitation_coeff.dot(pd.DataFrame(np.sqrt(self.t_electrons.value) *
                                                                       self.n_e.values).T)
        u0 = lines.nu.values[np.newaxis].T / self.t_electrons.value * (const.h.cgs.value / const.k_B.cgs.value)
        # mask_ions = lines.ion_number.values > 0
        # mask_ions = np.expand_dims(mask_ions, axis = 1) * np.ones(u0.shape[1], dtype=bool)
        # mask_neutral = np.logical_not(mask_ions)
        # TODO: gbar is 0.7 for transitions within a principal quantum number and is also different for neutral atoms
        # Cloudy uses a value of gbar = h * nu / (k_B * T_e)/10 for neutral atoms
        gamma = 0.276 * np.exp(u0) * expn(1, u0)
        # gamma[np.logical_and(gamma < 0.2, mask_ions)] = 0.2
        gamma[gamma < 0.2] = 0.2
        factor = pd.DataFrame(u0 * np.exp(-u0) * gamma)
        coll_excitation_coeff = coll_excitation_coeff.multiply(factor, axis=0)
        # coll_excitation_coeff.set_index(self.level_lower_index.set_names['atomic_number',
        # 'ion_number', 'source_level_number'], inplace=True)
        self._set_coll_excitation_coeff_index(coll_excitation_coeff, macro_atom_references)
        return coll_excitation_coeff

    def _set_coll_excitation_coeff_index(self, coll_excitation_coeff, macro_atom_references, ):
        source_level_idx = macro_atom_references.loc[self.level_lower_index, 'references_idx'].values
        destination_level_idx = macro_atom_references.loc[self.level_upper_index, 'references_idx'].values
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        coll_excitation_coeff.set_index(tmp_multi_index, inplace=True)

    def calculate_coll_deexcitation_coeff(self, coll_excitation_coeff, lte_level_pop):
        coll_deexcitation_coeff = coll_excitation_coeff.swaplevel(0, 1, axis=0)
        coll_deexcitation_coeff.index.names = coll_excitation_coeff.index.names[:]
        lte_level_pop_lower = lte_level_pop.loc[self.level_lower_index]
        lte_level_pop_upper = lte_level_pop.loc[self.level_upper_index]
        return coll_deexcitation_coeff * (lte_level_pop_lower.values / lte_level_pop_upper.values)

    def transition_probabilities(self, coll_excitation_coeff, coll_deexcitation_coeff, levels):
        # TODO: Check that c_einstein has the right units
        # TODO: Units seem to be off, check calculation of coll. coeff.
        c_einstein = (4. * (np.pi * const.e.gauss) ** 2 / (const.c.cgs * const.m_e.cgs)).value
        energy_lower = levels.loc[self.level_lower_index, 'energy'] * units.erg.to(units.eV)
        energy_upper = levels.loc[self.level_upper_index, 'energy'] * units.erg.to(units.eV)
        internal_jump_upwards = coll_excitation_coeff.multiply(energy_lower.values, axis=0) / c_einstein
        internal_jump_downwards = coll_deexcitation_coeff.multiply(energy_upper.values, axis=0) / c_einstein
        collisional_deexcitation = coll_deexcitation_coeff.multiply(
            (energy_upper.values - energy_lower.values), axis=0) / c_einstein
        return internal_jump_upwards, internal_jump_downwards, collisional_deexcitation

    def _calculate_excitation_cooling_rate(self, coll_excitation_coeff, levels, level_pop):
        energy_lower = levels.loc[self.level_lower_index, 'energy']
        energy_upper = levels.loc[self.level_upper_index, 'energy']
        cooling_rate = coll_excitation_coeff.multiply((energy_upper.values - energy_lower.values), axis=0)
        cooling_rate = cooling_rate.multiply(level_pop.loc[self.level_lower_index].values, axis=0)
        cooling_rate = cooling_rate.groupby(level=1).sum()
        return cooling_rate

    def _get_level_lower_index(self, lines):
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_lower']])

    def _get_level_upper_index(self, lines):
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_upper']])


class TransitionProbabilitiesCombined(object):
    def __init__(self, macro_atom_data, macro_atom_references, continuum_data,
                 radiative_prob, coll_int_jump_up_prob, coll_int_jump_down_prob, coll_deact_prob, ion_prob):
        self.combined_transition_probabilities = None
        self.block_references = None
        self.destination_level_id = None
        self.transition_probabilities_array = None
        self.radiative_probabilities = \
            self._prepare_radiative_probabilities(macro_atom_data=macro_atom_data, radiative_prob=radiative_prob,
                                                  macro_atom_references=macro_atom_references)
        self.collisional_deactivation_probabilities = \
            self._prepare_collisional_deactivation_probabilities(coll_deact_prob)
        self.ionization_probabilities = self._prepare_ionization_probabilities(
            ion_prob, macro_atom_references, continuum_data.continuum_references)
        # This is mostly useful when using Van Regemorter
        self.added_probabilities = self._add_internal_jump_probabilities(coll_int_jump_up_prob=coll_int_jump_up_prob,
                                                                         coll_int_jump_down_prob=coll_int_jump_down_prob,
                                                                         macro_atom_data=macro_atom_data)
        self._combine_transition_probabilities(self.added_probabilities, self.collisional_deactivation_probabilities,
                                               self.ionization_probabilities)

    def _prepare_radiative_probabilities(self, radiative_prob, macro_atom_data, macro_atom_references):
        source_level_idx = self._get_source_level_idx(macro_atom_data, macro_atom_references)
        destination_level_idx = macro_atom_data.destination_level_idx.values
        new_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                              names=['source_level_idx', 'destination_level_idx'])
        radiative_prob_prep = radiative_prob.set_index(new_index)
        radiative_prob_prep.insert(0, 'transition_type', macro_atom_data.transition_type.values)
        radiative_prob_prep.insert(1, 'lines_idx', macro_atom_data['lines_idx'].values)
        return radiative_prob_prep

    def _prepare_collisional_deactivation_probabilities(self, coll_deact_prob):
        coll_deact_prob_prep = coll_deact_prob.copy()
        coll_deact_prob_prep.insert(0, 'transition_type', -2 * np.ones(coll_deact_prob.values.shape[0]))
        coll_deact_prob_prep.insert(1, 'lines_idx', -1 * np.ones(coll_deact_prob.values.shape[0]))
        return coll_deact_prob_prep

    def _prepare_ionization_probabilities(self, ion_prob, macro_atom_references, macro_atom_continuum_references):
        # WARNING: destination level id is continuum id; the value itself is not unique
        ion_prob_prep = ion_prob.copy()
        ion_prob_prep.insert(0, 'transition_type', 2 * np.ones(ion_prob.values.shape[0]))
        ion_prob_prep.insert(1, 'lines_idx', -1 * np.ones(ion_prob.values.shape[0]))
        multi_index = self._get_ion_prob_index(ion_prob_prep, macro_atom_references, macro_atom_continuum_references)
        ion_prob_prep.set_index(multi_index, inplace=True)
        return ion_prob_prep

    def _get_source_level_idx(self, macro_atom_data, macro_atom_references):
        source_level_index = pd.MultiIndex.from_arrays([macro_atom_data['atomic_number'], macro_atom_data['ion_number'],
                                                        macro_atom_data['source_level_number']])
        return macro_atom_references.loc[source_level_index, 'references_idx'].values

    def _add_internal_jump_probabilities(self, coll_int_jump_up_prob, coll_int_jump_down_prob, macro_atom_data):
        added_probabilities = self.radiative_probabilities.copy()
        transition_up_filter = (macro_atom_data.transition_type == 1).values
        added_probabilities[transition_up_filter] = \
            added_probabilities[transition_up_filter].combineAdd(coll_int_jump_up_prob)
        transition_down_filter = (macro_atom_data.transition_type == 0).values
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
        self.transition_probabilities_array = np.ascontiguousarray(
            self.combined_transition_probabilities.ix[:, 2:].values)
        self.transition_type = self.combined_transition_probabilities['transition_type'].values.astype(np.int64)
        self.transition_line_id = self.combined_transition_probabilities['lines_idx'].values.astype(np.int64)

    def _get_new_block_references(self, combined_probabilities):
        block_references = combined_probabilities[0].groupby(level=0).count().cumsum().values
        block_references = np.hstack([[0], block_references])
        return block_references

    def _get_ion_prob_index(self, ion_prob, macro_atom_references, macro_atom_continuum_references):
        level_lower_index = ion_prob.index
        source_level_idx = macro_atom_references.loc[level_lower_index, 'references_idx'].values
        destination_level_idx = self._get_continuum_idx(level_lower_index, macro_atom_continuum_references)
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        return tmp_multi_index

    def _get_continuum_idx(self, multi_index_full, continuum_references):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1)
        ion_number_index = pd.MultiIndex.from_arrays([atomic_number, ion_number])
        return continuum_references.loc[ion_number_index, 'references_idx']


class CoolingRates(object):
    def __init__(self, plasma_array, coll_excitation_cooling_rate, coll_ionization_cooling_rate, fb_cooling_rate):
        self.ff_cooling = \
            self._calculate_ff_cooling_rate(electron_densities=plasma_array.electron_densities,
                                            t_electrons=plasma_array.t_electrons,
                                            ion_number_density=plasma_array.ion_number_density,
                                            ion_charges=plasma_array.ion_number_density.index.get_level_values(
                                                1).values)

        self.collisional_excitation_cooling = coll_excitation_cooling_rate
        self.collisional_ionization_cooling = coll_ionization_cooling_rate
        self.fb_cooling = fb_cooling_rate
        self._set_total_cooling_rates()
        self._set_cooling_probabilities()
        self._set_individual_cooling_probabilities()
        self._prepare_montecarlo_data()

    def _calculate_ff_cooling_rate(self, electron_densities, t_electrons, ion_number_density, ion_charges):
        C0 = 1.426e-27  # in cgs units (see Osterbrock 1974, newer version C0 = 1.42)
        # TODO: check that units are used consistently
        # TODO: value for Gaunt factor (Lucy: = 1; Osterbrock recommendation for nebular conditions: = 1.3 )
        factor = ion_number_density.mul(np.square(ion_charges), axis=0).sum().values
        cooling_rate = C0 * electron_densities * np.sqrt(t_electrons) * factor
        return cooling_rate.values

    def _set_total_cooling_rates(self):
        self.collisional_excitation_cooling_total = self.collisional_excitation_cooling.sum().values
        self.collisional_ionization_cooling_total = self.collisional_ionization_cooling.sum().values
        self.fb_cooling_total = self.fb_cooling.sum().values

    def _set_cooling_probabilities(self):
        total_cooling_rate = self.fb_cooling_total + self.ff_cooling + self.collisional_excitation_cooling_total + \
                             self.collisional_ionization_cooling_total
        self.fb_cooling_prob = self.fb_cooling_total / total_cooling_rate
        self.ff_cooling_prob = self.ff_cooling / total_cooling_rate
        self.coll_ion_cooling_prob = self.collisional_ionization_cooling_total / total_cooling_rate
        self.coll_exc_cooling_prob = self.collisional_excitation_cooling_total / total_cooling_rate

    def _set_individual_cooling_probabilities(self):
        self.fb_cooling_prob_individual = self.fb_cooling.divide(self.fb_cooling_total, axis=1)
        self.coll_exc_cooling_prob_individual = self.collisional_excitation_cooling.divide(
            self.collisional_excitation_cooling_total, axis=1)
        self.coll_ion_cooling_prob_individual = self.collisional_ionization_cooling.divide(
            self.collisional_ionization_cooling_total, axis=1)

    def _prepare_montecarlo_data(self):
        self.fb_cooling_prob_array = self._get_contiguous_array(self.fb_cooling_prob_individual)
        self.coll_exc_cooling_prob_array = self._get_contiguous_array(self.coll_exc_cooling_prob_individual)
        self.coll_ion_cooling_prob_array = self._get_contiguous_array(self.coll_ion_cooling_prob_individual)

    def _get_contiguous_array(self, dataframe):
        return np.ascontiguousarray(dataframe.values.transpose())

    @property
    def fb_cooling_prob_nd(self):
        return self.fb_cooling_prob_array.shape[1]

    @property
    def coll_exc_cooling_prob_nd(self):
        return self.coll_exc_cooling_prob_array.shape[1]

    @property
    def coll_ion_cooling_prob_nd(self):
        return self.coll_ion_cooling_prob_array.shape[1]





