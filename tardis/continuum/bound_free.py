import numpy as np
import logging
from astropy import units
from astropy import constants as const
import pandas as pd
from tardis.util import intensity_black_body
from scipy.integrate import simps
from tardis import macro_atom
from tardis.continuum.base import ContinuumProcess

logger = logging.getLogger(__name__)


class IonizationRates(ContinuumProcess):
    def __init__(self, input_data):
        super(IonizationRates, self).__init__(input_data)

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
        j_nus = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads)
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
        collion_coeff = collion_coeff.multiply(self.electron_densities / np.sqrt(self.t_electrons), axis=1)
        return collion_coeff

    def _calculate_boltzmann_factor(self, nu):
        u0s = self._calculate_u0s(nu)
        return np.exp(-u0s)

    def _calculate_u0s(self, nu):
        u0s = nu[np.newaxis].T / self.t_electrons * (const.h.cgs.value / const.k_B.cgs.value)
        return u0s

    def _calculate_factor(self, nu_ijk, index):
        u0s = self._calculate_u0s(nu_ijk)
        factor = 1. / u0s * np.exp(-u0s)
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
        source_level_index = self.collisional_ionization_coefficient.index
        level_pop = self._get_level_pop(source_level_index)
        coll_ion_cooling_rate = self.collisional_ionization_coefficient.multiply(level_pop)
        coll_ion_cooling_rate = coll_ion_cooling_rate.multiply(self.nu_i * const.h.cgs.value, axis=0)
        continuum_idx = self._get_continuum_idx(coll_ion_cooling_rate.index)
        coll_ion_cooling_rate.set_index(continuum_idx, inplace=True)
        return coll_ion_cooling_rate

    @property
    def transition_probabilities2cont(self):
        trans_prob2cont = \
            self.corrected_photoionization_coefficient.add(self.collisional_ionization_coefficient) / self.c_einstein
        # trans_prob2cont = trans_prob2cont.multiply(self.nu_i, axis=0) * const.h.cgs.value * units.erg.to(units.eV)
        source_level_energy = self._get_level_energy(trans_prob2cont.index) * units.erg.to(units.eV)
        trans_prob2cont = trans_prob2cont.multiply(source_level_energy, axis=0)
        return trans_prob2cont


class RadiativeRecombination(ContinuumProcess):
    def __init__(self, input_data):
        super(RadiativeRecombination, self).__init__(input_data)

        self.block_references = None
        self.sp_recombination_coeff = self._calculate_sp_recombination_coeff()
        self.transition_probabilities = self._calculate_transition_probabilities()
        self.fb_cooling_rate = self._calculate_fb_cooling_rate()
        self._set_montecarlo_data()

    def _calculate_sp_recombination_coeff(self, modified=False):
        if modified == False:
            recomb_coeff = (8 * np.pi * self.photoionization_data['x_sect']
                            * (self.photoionization_data['nu']) ** 2 / (const.c.cgs.value) ** 2).values
        else:
            recomb_coeff = (8 * np.pi * self.photoionization_data['x_sect']
                            * (self.photoionization_data['nu']) ** 3 / (const.c.cgs.value) ** 2).values

        recomb_coeff = recomb_coeff[:, np.newaxis]
        boltzmann_factor = np.exp(-self.photoionization_data.nu.values[np.newaxis].T / \
                                  self.t_rads * (const.h.cgs.value / const.k_B.cgs.value))
        recomb_coeff = pd.DataFrame(boltzmann_factor * recomb_coeff, index=self.photoionization_data.index)
        recomb_coeff = recomb_coeff.divide(self.electron_densities, axis=1)
        recomb_coeff.insert(0, 'nu', self.photoionization_data['nu'])
        recomb_coeff = recomb_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(self.no_of_shells):
            tmp[i] = recomb_coeff.apply(lambda sub: simps(sub[i], sub['nu'], even='first'))
            if modified == True:
                tmp[i] /= self.nu_i
        recomb_coeff = pd.DataFrame(tmp)
        recomb_coeff = recomb_coeff.multiply(self._get_lte_level_pop(recomb_coeff.index))
        ion_number_density = self._get_ion_number_density(recomb_coeff.index)
        recomb_coeff = recomb_coeff.divide(ion_number_density.values)
        return recomb_coeff

    def _calculate_transition_probabilities(self):
        internal_jump_prob = self._calculate_internal_jump_probabilities()
        deactivation_prob = self._calculate_deactivation_probabilities()
        trans_prob = pd.concat([internal_jump_prob, deactivation_prob])

        block_references = self._get_block_references(trans_prob)
        self.block_references = block_references
        trans_prob = self._normalize_transition_probabilities(trans_prob, no_ref_columns=3)
        return trans_prob

    def _calculate_internal_jump_probabilities(self):
        target_level_energy = self._get_level_energy(self.sp_recombination_coeff.index)
        trans_prob = self.sp_recombination_coeff.multiply(target_level_energy, axis=0)

        destination_level_idx = self._get_level_idx(trans_prob.index)
        trans_prob.insert(0, 'destination_level_idx', destination_level_idx)
        trans_prob.insert(1, 'continuum_edge_idx', - 1 * np.ones(trans_prob.shape[0], dtype=np.int64))
        trans_prob.insert(2, 'transition_type', np.zeros(trans_prob.shape[0], dtype=np.int64))
        return trans_prob

    def _calculate_deactivation_probabilities(self):
        energy_difference = self.nu_i * const.h.cgs.value
        trans_prob = self.sp_recombination_coeff.multiply(energy_difference, axis=0)

        continuum_edge_idx = self._get_continuum_edge_idx(trans_prob.index)
        trans_prob.insert(0, 'destination_level_idx', - 1 * np.ones(trans_prob.shape[0], dtype=np.int64))
        trans_prob.insert(1, 'continuum_edge_idx', continuum_edge_idx)
        trans_prob.insert(2, 'transition_type', -3 * np.ones(trans_prob.shape[0], dtype=np.int64))
        return trans_prob

    def _get_ion_number_density(self, multi_index_full):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1) + 1
        ion_number_index = pd.MultiIndex.from_arrays([atomic_number, ion_number])
        return self.ion_number_density.loc[ion_number_index]

    def _calculate_fb_cooling_rate(self):
        self.sp_recombination_coeff_E = self._calculate_sp_recombination_coeff(modified=True)
        fb_cooling_rate = (self.sp_recombination_coeff_E - self.sp_recombination_coeff)
        fb_cooling_rate = fb_cooling_rate.multiply(const.h.cgs.value * self.nu_i, axis=0)
        fb_cooling_rate = fb_cooling_rate.multiply(self.electron_densities, axis=1)
        ion_number_density = self._get_ion_number_density(fb_cooling_rate.index)
        fb_cooling_rate = fb_cooling_rate.multiply(ion_number_density.values)
        continuum_edge_idx = self._get_continuum_edge_idx(fb_cooling_rate.index)
        fb_cooling_rate.set_index(continuum_edge_idx, inplace=True)
        return fb_cooling_rate

    def _set_montecarlo_data(self):
        self.data_array = self._get_contiguous_array(self.transition_probabilities.ix[:, 3:])
        self.data_array_nd = self.data_array.shape[1]

    # To be depreciated
    @property
    def data(self):
        return self.transition_probabilities

    def _calculate_transition_probabilities(self):
        import ipdb;

        ipdb.set_trace()
        trans_prob = pd.concat([self.sp_recombination_coeff, self.sp_recombination_coeff])
        # TODO: In the future, we should check if the photoionization_data and the macro_atom_continuum_data have the
        # same structure (maybe do this in the preparation of the continuum_data)
        trans_prob = trans_prob.multiply(self.input.macro_atom_continuum_data.transition_probability.values
                                         * units.eV.to(units.erg), axis=0)
        # WARNING: Not sure if this is safe under all circumstances
        block_references = self._get_block_references(trans_prob)
        macro_atom.normalize_transition_probabilities(trans_prob.values, block_references)
        trans_prob.insert(0, 'destination_level_idx',
                          self.input.macro_atom_continuum_data.level_lower_idx.values)
        trans_prob.insert(1, 'continuum_edge_idx',
                          self.input.macro_atom_continuum_data.continuum_edge_idx.values)
        trans_prob.insert(2, 'transition_type',
                          self.input.macro_atom_continuum_data.transition_type.values)
        return trans_prob