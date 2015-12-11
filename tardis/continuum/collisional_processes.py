import numpy as np
import pandas as pd
import astropy.constants as const
from scipy.special import expn

from tardis.continuum.base import PhysicalContinuumProcess, InverseProcess, BoundFreeEnergyMixIn
from tardis.continuum.constants import continuum_constants as cconst


class CollisionalExcitation(PhysicalContinuumProcess):
    name = 'collisional_excitation'
    macro_atom_transitions = 'up'

    def __init__(self, input_data, mode='Van Regemorter'):
        self.mode = mode
        super(CollisionalExcitation, self).__init__(input_data)

    def _calculate_rate_coefficient(self):
        if self.mode == 'Van Regemorter':
            return self._calculate_rate_coefficient_regemorter()
        else:
            raise NotImplementedError

    def _calculate_cooling_rate(self):
        energy_lower = self._get_level_energy(self.level_lower_index)
        energy_upper = self._get_level_energy(self.level_upper_index)
        level_pop_lower = self._get_level_pop(self.level_lower_index)

        cooling_rate = self.rate_coefficient.multiply(energy_upper - energy_lower, axis=0)
        cooling_rate = cooling_rate.multiply(level_pop_lower, axis=0)
        cooling_rate = cooling_rate.groupby(level=1).sum()
        return cooling_rate

    def _calculate_rate_coefficient_regemorter(self):
        self.level_lower_index = self._get_level_lower_index()
        self.level_upper_index = self._get_level_upper_index()
        I_H = cconst.I_H
        f_lu = self.input.lines.f_lu.values
        n_e = self.electron_densities
        nu_lines = self.input.lines.nu.values

        coll_excitation_coeff = pd.DataFrame(14.5 * cconst.c0_reg * f_lu * (I_H / (const.h.cgs.value * nu_lines)) ** 2)
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

    def _get_level_lower_index(self):
        lines = self.input.lines
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_lower']])

    def _get_level_upper_index(self):
        lines = self.input.lines
        return pd.MultiIndex.from_arrays([lines['atomic_number'], lines['ion_number'], lines['level_number_upper']])

    @property
    def level_lower_energy(self):
        return self._get_level_energy(self.level_lower_index)

    @property
    def level_upper_energy(self):
        return self._get_level_energy(self.level_upper_index)


class CollisionalDeexcitation(InverseProcess):
    name = 'collisional_deexcitation'
    name_of_inverse_process = 'collisional_excitation'
    macro_atom_transitions = 'down'

    @classmethod
    def _calculate_inverse_rate(cls, inverse_process):
        rate_coefficient = inverse_process.rate_coefficient.swaplevel(0, 1, axis=0)
        rate_coefficient.index.names = inverse_process.rate_coefficient.index.names[:]
        lte_level_pop_lower = inverse_process._get_lte_level_pop(inverse_process.level_lower_index)
        lte_level_pop_upper = inverse_process._get_lte_level_pop(inverse_process.level_upper_index)
        return rate_coefficient * (lte_level_pop_lower / lte_level_pop_upper)


class CollisionalRecombination(InverseProcess):
    name = 'collisional_recombination'
    name_of_inverse_process = 'collisional_ionization'

    @classmethod
    def _calculate_inverse_rate(cls, inverse_process):
        level_lower_index = inverse_process.rate_coefficient.index
        lte_level_pop_lower = inverse_process._get_lte_level_pop(level_lower_index)
        lte_pop_continuum = inverse_process._get_ion_number_density(level_lower_index)
        rate_coeff = inverse_process.rate_coefficient * (lte_level_pop_lower / lte_pop_continuum)
        # Needed for consistency with the radiative recombination rate coefficient
        rate_coeff = rate_coeff.divide(inverse_process.electron_densities, axis=1)
        return rate_coeff


class CollisionalIonization(PhysicalContinuumProcess, BoundFreeEnergyMixIn):
    name = 'collisional_ionization'
    macro_atom_transitions = 'continuum'

    def __init__(self, input_data):
        super(CollisionalIonization, self).__init__(input_data)
        # self._set_ionization_rates_index(self.rate_coefficient)

    def _calculate_rate_coefficient(self):
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

    def _calculate_cooling_rate(self):
        source_level_index = self.rate_coefficient.index
        level_pop = self._get_level_pop(source_level_index)
        coll_ion_cooling_rate = self.rate_coefficient.multiply(level_pop)
        coll_ion_cooling_rate = coll_ion_cooling_rate.multiply(self.nu_i * const.h.cgs.value, axis=0)
        continuum_idx = self._get_continuum_idx(coll_ion_cooling_rate.index)
        coll_ion_cooling_rate.set_index(continuum_idx, inplace=True)
        coll_ion_cooling_rate.index.set_names('continuum_idx')
        return coll_ion_cooling_rate

    def _calculate_factor(self, nu_ijk, index):
        u0s = self._calculate_u0s(nu_ijk)
        factor = 1. / u0s * np.exp(-u0s)
        factor = pd.DataFrame(factor, index=index, columns=np.arange(len(self.t_electrons)))
        return factor