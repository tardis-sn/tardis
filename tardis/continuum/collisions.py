from tardis.continuum.base import ContinuumProcess
import numpy as np
import pandas as pd
import astropy.constants as const
import astropy.units as units
from scipy.special import expn


class CollisionalExcitation(ContinuumProcess):
    def __init__(self, input_data, mode='Van Regemorter'):
        self.input = input_data
        self.mode = mode
        if mode == 'Van Regemorter':
            self.level_lower_index = self._get_level_lower_index()
            self.level_upper_index = self._get_level_upper_index()
            self.coll_excitation_coeff = self.calculate_coll_excitation_coeff_regemorter()
        else:
            raise NotImplementedError
        self.excitation_cooling_rate = self._calculate_excitation_cooling_rate(self.coll_excitation_coeff)

        self.coll_deexcitation_coeff = self.calculate_coll_deexcitation_coeff()
        self.transition_prob_int_up, self.transition_prob_int_down, self.transition_prob_deact = \
            self.transition_probabilities(coll_excitation_coeff=self.coll_excitation_coeff,
                                          coll_deexcitation_coeff=self.coll_deexcitation_coeff)

    def calculate_coll_excitation_coeff_regemorter(self):
        I_H = self.input.I_H
        f_lu = self.input.lines.f_lu.values
        n_e = self.electron_densities
        nu_lines = self.input.lines.nu.values

        coll_excitation_coeff = pd.DataFrame(14.5 * self.c0_reg * f_lu * (I_H / (const.h.cgs.value * nu_lines)) ** 2)
        coll_excitation_coeff = coll_excitation_coeff.dot(pd.DataFrame(np.sqrt(self.t_electrons) * n_e).T)
        u0 = nu_lines[np.newaxis].T / self.t_electrons * (const.h.cgs.value / const.k_B.cgs.value)
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
        self._set_coll_excitation_coeff_index(coll_excitation_coeff)
        return coll_excitation_coeff

    def _set_coll_excitation_coeff_index(self, coll_excitation_coeff):
        source_level_idx = self._get_level_idx(self.level_lower_index)
        destination_level_idx = self._get_level_idx(self.level_upper_index)
        tmp_multi_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                                    names=['source_level_idx', 'destination_level_idx'])
        coll_excitation_coeff.set_index(tmp_multi_index, inplace=True)

    def calculate_coll_deexcitation_coeff(self):
        coll_deexcitation_coeff = self.coll_excitation_coeff.swaplevel(0, 1, axis=0)
        coll_deexcitation_coeff.index.names = self.coll_excitation_coeff.index.names[:]
        lte_level_pop_lower = self._get_lte_level_pop(self.level_lower_index)
        lte_level_pop_upper = self._get_lte_level_pop(self.level_upper_index)
        return coll_deexcitation_coeff * (lte_level_pop_lower / lte_level_pop_upper)

    def transition_probabilities(self, coll_excitation_coeff, coll_deexcitation_coeff):
        energy_lower = self._get_level_energy(self.level_lower_index) * units.erg.to(units.eV)
        energy_upper = self._get_level_energy(self.level_upper_index) * units.erg.to(units.eV)

        internal_jump_upwards = coll_excitation_coeff.multiply(energy_lower, axis=0) / self.c_einstein
        internal_jump_downwards = coll_deexcitation_coeff.multiply(energy_upper, axis=0) / self.c_einstein
        collisional_deexcitation = coll_deexcitation_coeff.multiply(
            (energy_upper - energy_lower), axis=0) / self.c_einstein
        return internal_jump_upwards, internal_jump_downwards, collisional_deexcitation

    def _calculate_excitation_cooling_rate(self, coll_excitation_coeff):
        energy_lower = self._get_level_energy(self.level_lower_index)
        energy_upper = self._get_level_energy(self.level_upper_index)
        level_pop_lower = self._get_level_pop(self.level_lower_index)

        cooling_rate = coll_excitation_coeff.multiply(energy_upper - energy_lower, axis=0)
        cooling_rate = cooling_rate.multiply(level_pop_lower, axis=0)
        cooling_rate = cooling_rate.groupby(level=1).sum()
        return cooling_rate

    def _get_level_lower_index(self):
        lines = self.input.lines
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_lower']])

    def _get_level_upper_index(self):
        lines = self.input.lines
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_upper']])