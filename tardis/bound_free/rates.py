import numpy as np
import logging
import os
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
    def __init__(self, photoionization_data, ws, t_rads, electron_densities):
        self.photoionization_data = photoionization_data.copy()
        self.ws = ws
        self.t_rads = t_rads
        self.electron_densities = electron_densities
        self.corrected_photoionization_coefficient = self._calculate_corrected_photoion_coeff()

    # TODO: Stimulated emission is missing
    def _calculate_corrected_photoion_coeff(self):
        j_nus = self._calculate_j_nus()
        corrected_photoion_coeff = j_nus.multiply(4. * np.pi * self.photoionization_data['x_sect'] /
                                                  self.photoionization_data['nu'] / const.h.cgs.value, axis=0)
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

    # TODO: This needs some refactoring
    def calculate_collion_coeff(self):
        collion_coeff = self.photoionization_data.groupby(level=[0, 1, 2]).first()
        factor = self._calculate_factor(collion_coeff['nu'].values)
        factor = pd.DataFrame(factor, index=collion_coeff.index, columns=np.arange(len(self.t_rads)))
        collion_coeff = 1.55e13 * collion_coeff['x_sect']
        # TODO: Try doing this without resetting the index
        collion_coeff = collion_coeff.reset_index()
        collion_coeff.ix[collion_coeff['ion_number'] == 0, 'x_sect'] *= 0.1
        collion_coeff.ix[collion_coeff['ion_number'] == 1, 'x_sect'] *= 0.2
        collion_coeff.ix[collion_coeff['ion_number'] >= 2, 'x_sect'] *= 0.3
        collion_coeff.set_index(['atomic_number', 'ion_number', 'level_number'], inplace=True)
        collion_coeff = factor.multiply(collion_coeff.x_sect, axis=0)
        collion_coeff = collion_coeff.multiply(self.electron_densities / np.sqrt(0.9 * self.t_rads.value), axis=1)
        return collion_coeff

    def _calculate_factor(self, nu_ijk):
        # WARNING: This might be unsafe for future implementations.
        t_electrons = 0.9 * self.t_rads
        u0s = nu_ijk[np.newaxis].T / t_electrons.value * (const.h.cgs.value / const.k_B.cgs.value)
        factor = 1 / u0s * np.exp(-u0s)
        return factor

    # TODO: Check that the right units are used
    @property
    def transition_probabilities2cont(self):
        # TODO: multiply with energy difference
        c_einstein = (4. * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)).value
        return self.corrected_photoionization_coefficient / c_einstein


class CollisionalRates(object):
    def __init__(self, lines, t_electrons, electron_densities, lte_level_pop,
                 levels, macro_atom_references, mode='Van Regemorter'):
        # TODO: t_electrons and n_e should not be attributes
        self.t_electrons = t_electrons
        self.n_e = electron_densities
        self.mode = mode
        if mode == 'Van Regemorter':
            self.level_lower_index = pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'],
                                                                lines['level_number_lower']])
            self.level_upper_index = pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'],
                                                                lines['level_number_upper']])
            self.coll_excitation_coeff = self.calculate_coll_excitation_coeff_regemorter(lines=lines,
                                                                                         macro_atom_references=macro_atom_references)
        else:
            raise NotImplementedError

        self.collisional_bb_cooling = self.bb_cooling_rate(level_pop=lte_level_pop, levels=levels,
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
        source_level_idx = macro_atom_references.loc[self.level_lower_index, 'references_idx'].values
        destination_level_idx = macro_atom_references.loc[self.level_upper_index, 'references_idx'].values
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        coll_excitation_coeff.set_index(tmp_multi_index, inplace=True)
        return coll_excitation_coeff

    def calculate_coll_deexcitation_coeff(self, coll_excitation_coeff, lte_level_pop):
        coll_deexcitation_coeff = coll_excitation_coeff.swaplevel(0, 1, axis=0)
        coll_deexcitation_coeff.index.rename(coll_excitation_coeff.index.names, inplace=True)
        level_pop_lower = lte_level_pop.loc[self.level_lower_index]
        level_pop_upper = lte_level_pop.loc[self.level_upper_index]
        return coll_excitation_coeff * (level_pop_lower.values / level_pop_upper.values)

    def transition_probabilities(self, coll_excitation_coeff, coll_deexcitation_coeff, levels):
        # TODO: Check that c_einstein has the right units
        # TODO: Units seem to be off, check calculation of coll. coeff.
        c_einstein = (4. * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)).value
        energy_lower = levels.loc[self.level_lower_index, 'energy'] * units.erg.to(units.eV)
        energy_upper = levels.loc[self.level_upper_index, 'energy'] * units.erg.to(units.eV)
        internal_jump_upwards = coll_excitation_coeff.multiply(energy_lower.values, axis=0) / c_einstein
        internal_jump_downwards = coll_deexcitation_coeff.multiply(energy_upper.values, axis=0) / c_einstein
        collisional_deexcitation = coll_deexcitation_coeff.multiply(
            (energy_upper.values - energy_lower.values), axis=0) / c_einstein
        return internal_jump_upwards, internal_jump_downwards, collisional_deexcitation

    def bb_cooling_rate(self, coll_excitation_coeff, levels, level_pop):
        # TODO: use consistent units
        energy_lower = levels.loc[self.level_lower_index, 'energy']
        energy_upper = levels.loc[self.level_upper_index, 'energy']
        cooling_rate = coll_excitation_coeff.multiply((energy_upper.values - energy_lower.values), axis=0)
        cooling_rate = cooling_rate.multiply(level_pop.loc[self.level_lower_index].values, axis=0)
        cooling_rate = cooling_rate.groupby(level=1).sum()
        return cooling_rate


class TransitionProbabilitiesCombined(object):
    def __init__(self, macro_atom_data, macro_atom_references, radiative_prob, coll_excitation_prob):
        self.radiative_probabilities = \
            self._prepare_radiative_probabilities(macro_atom_data=macro_atom_data, radiative_prob=radiative_prob,
                                                  macro_atom_references=macro_atom_references)
        # This is mostly useful when using Van Regemorter
        # self._add_internal_jump_probabilities(coll_excitation_prob=coll_excitation_prob,
        #                                      macro_atom_data=macro_atom_data)

    def _prepare_radiative_probabilities(self, radiative_prob, macro_atom_data, macro_atom_references):
        source_level_idx = self._get_source_level_idx(macro_atom_data, macro_atom_references)
        destination_level_idx = macro_atom_data.destination_level_idx.values
        new_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                              names=[u'source_level_idx', u'destination_level_idx'])
        radiative_prob_prep = radiative_prob.set_index(new_index)
        # TODO: insert line_idx and transition_type
        return radiative_prob_prep

    def _get_source_level_idx(self, macro_atom_data, macro_atom_references):
        source_level_index = pd.MultiIndex.from_arrays([macro_atom_data['atomic_number'], macro_atom_data['ion_number'],
                                                        macro_atom_data['source_level_number']])
        return macro_atom_references.loc[source_level_index, 'references_idx'].values

    def _add_internal_jump_probabilities(self, coll_excitation_prob, macro_atom_data):
        added_probabilities = self.radiative_probabilities.copy()
        transition_up_filter = (macro_atom_data.transition_type == 1).values
        added_probabilities[transition_up_filter] = \
            added_probabilities[transition_up_filter].add(coll_excitation_prob)
        return added_probabilities

    def _get_new_block_references(self, combined_probabilities):
        block_references = combined_probabilities[0].groupby(level=0).count().cumsum().values
        block_references = np.hstack([[0], block_references])
        return block_references






