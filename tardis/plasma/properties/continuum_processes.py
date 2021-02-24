import logging

import numpy as np
import pandas as pd

from numba import prange, njit
from astropy import constants as const

from tardis.plasma.exceptions import PlasmaException
from tardis.plasma.properties.base import (ProcessingPlasmaProperty,
                                           Input,
                                           TransitionProbabilitiesProperty)
from tardis.plasma.properties.j_blues import JBluesDiluteBlackBody

__all__ = ['SpontRecombRateCoeff', 'StimRecombRateCoeff', 'PhotoIonRateCoeff',
           'PhotoIonEstimatorsNormFactor', 'PhotoIonRateCoeffEstimator',
           'StimRecombRateCoeffEstimator', 'CorrPhotoIonRateCoeff',
           'BfHeatingRateCoeffEstimator', 'SpontRecombCoolingRateCoeff',
           'BaseRecombTransProbs', 'BasePhotoIonTransProbs',
           'CollDeexcRateCoeff', 'CollExcRateCoeff', 'BaseCollisionTransProbs',
           'AdiabaticCoolingRate', 'FreeFreeCoolingRate',
           'BoundFreeOpacity', 'LevelNumberDensityLTE',
           'PhotoIonBoltzmannFactor']

logger = logging.getLogger(__name__)

njit_dict = {'fastmath': False, 'parallel': False}


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrate a function over blocks.

    This function integrates a function `f` defined at locations `x`
    over blocks given in `block_references`.

    Parameters
    ----------
    f : Two-dimensional Numpy Array, dtype float
        Input array to integrate.
    x : One-dimensional Numpy Array, dtype float
        The sample points corresponding to the `f` values.
    block_references : One-dimensional Numpy Array, dtype int
        The start indices of the blocks to be integrated.

    Returns
    -------
    Two-dimensional Numpy Array, dtype float
        Array with integrated values.
    """
    integrated = np.zeros((len(block_references) - 1, f.shape[1]))
    for i in prange(f.shape[1]):  # columns
        for j in prange(len(integrated)):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[j, i] = np.trapz(f[start:stop, i], x[start:stop])
    return integrated


def get_ion_multi_index(multi_index_full, next_higher=True):
    """
    Calculate the corresponding ion MultiIndex for a level MultiIndex.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)
    next_higher : bool, default True
        If True use ion number of next higher ion, else use ion_number from
        multi_index_full.

    Returns
    -------
    Pandas MultiIndex (atomic_number, ion_number)
       Ion MultiIndex for the given level MultiIndex.
    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


def get_ground_state_multi_index(multi_index_full):
    """
    Calculate the ground-state MultiIndex for the next higher ion.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)

    Returns
    -------
    Pandas MultiIndex (atomic_number, ion_number)
        Ground-state MultiIndex for the next higher ion.
    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1) + 1
    level_number = np.zeros_like(ion_number)
    return pd.MultiIndex.from_arrays([atomic_number, ion_number, level_number])


def cooling_rate_series2dataframe(cooling_rate_series,
                                  destination_level_idx):
    """
    Transform cooling-rate Series to DataFrame.

    This function transforms a Series with cooling rates into
    an indexed DataFrame that can be used in MarkovChainTransProbs.

    Parameters
    ----------
    cooling_rate_series : Pandas Series, dtype float
        Cooling rates for a process with a single destination idx.
        Examples are adiabatic cooling or free-free cooling.
    destination_level_idx: String
        Destination idx of the cooling process; for example
        'adiabatic' for adiabatic cooling.

    Returns
    -------
    cooling_rate_frame : Pandas DataFrame, dtype float
        Indexed by source_level_idx, destination_level_idx, transition_type
        for the use in MarkovChainTransProbs.
    """
    index_names = ['source_level_idx', 'destination_level_idx',
                   'transition_type']
    index = pd.MultiIndex.from_tuples(
        [('k', destination_level_idx, -1)], names=index_names
    )
    cooling_rate_frame = pd.DataFrame(
        cooling_rate_series.values[np.newaxis],
        index=index
    )
    return cooling_rate_frame


class IndexSetterMixin(object):
    @staticmethod
    def set_index(p, photo_ion_idx, transition_type=0, reverse=True):
        idx = photo_ion_idx.loc[p.index]
        transition_type = transition_type * np.ones_like(
            idx.destination_level_idx
        )
        transition_type = pd.Series(transition_type, name='transition_type')
        idx_arrays = [idx.source_level_idx, idx.destination_level_idx]
        if reverse:
            idx_arrays = idx_arrays[::-1]
        idx_arrays.append(transition_type)
        index = pd.MultiIndex.from_arrays(idx_arrays)
        if reverse:
            index.names = index.names[:-1][::-1] + [index.names[-1]]
        p = p.set_index(index, drop=True)
        return p


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : Pandas DataFrame, dtype float
               The rate coefficient for spontaneous recombination.
    """
    outputs = ('alpha_sp',)
    latex_name = ('\\alpha^{\\textrm{sp}}',)

    def calculate(self, photo_ion_cross_sections, t_electrons,
                  photo_ion_block_references, photo_ion_index, phi_ik,
                  boltzmann_factor_photo_ion):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values

        alpha_sp = (8 * np.pi * x_sect * nu ** 2 / (const.c.cgs.value) ** 2)
        alpha_sp = alpha_sp[:, np.newaxis]
        alpha_sp = alpha_sp * boltzmann_factor_photo_ion
        alpha_sp = integrate_array_by_blocks(alpha_sp, nu,
                                             photo_ion_block_references)
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class SpontRecombCoolingRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    c_fb_sp : Pandas DataFrame, dtype float
              The rate coefficient for cooling by
              spontaneous recombination.
    """
    outputs = ('c_fb_sp',)
    latex_name = ('ca^{\\textrm{sp}}_{\\textrm{fb}}',)

    def calculate(self, photo_ion_cross_sections, t_electrons,
                  photo_ion_block_references, photo_ion_index, phi_ik, nu_i,
                  boltzmann_factor_photo_ion):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values
        factor = (1 - nu_i / photo_ion_cross_sections['nu']).values
        alpha_sp = (8 * np.pi * x_sect * factor * nu ** 3 /
                    (const.c.cgs.value) ** 2) * const.h.cgs.value
        alpha_sp = alpha_sp[:, np.newaxis]
        alpha_sp = alpha_sp * boltzmann_factor_photo_ion
        alpha_sp = integrate_array_by_blocks(alpha_sp, nu,
                                             photo_ion_block_references)
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
            The rate coefficient for radiative ionization.
    """
    outputs = ('gamma',)
    latex_name = ('\\gamma',)

    def calculate(self, photo_ion_cross_sections, gamma_estimator,
                  photo_ion_norm_factor, photo_ion_block_references,
                  photo_ion_index, t_rad, w):
        # Used for initialization
        if gamma_estimator is None:
            gamma = self.calculate_from_dilute_bb(photo_ion_cross_sections,
                                                  photo_ion_block_references,
                                                  photo_ion_index, t_rad, w)
        else:
            gamma = gamma_estimator * photo_ion_norm_factor
        return gamma

    @staticmethod
    def calculate_from_dilute_bb(photo_ion_cross_sections,
                                 photo_ion_block_references,
                                 photo_ion_index, t_rad, w):
        nu = photo_ion_cross_sections['nu']
        x_sect = photo_ion_cross_sections['x_sect']
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        gamma = j_nus.multiply(
            4. * np.pi * x_sect / nu / const.h.cgs.value, axis=0
        )
        gamma = integrate_array_by_blocks(gamma.values, nu.values,
                                          photo_ion_block_references)
        gamma = pd.DataFrame(gamma, index=photo_ion_index)
        return gamma


class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_stim : Pandas DataFrame, dtype float
                 The rate coefficient for stimulated recombination.
    """
    outputs = ('alpha_stim',)
    latex_name = ('\\alpha^{\\textrm{stim}}',)

    def calculate(self, photo_ion_cross_sections, alpha_stim_estimator,
                  photo_ion_norm_factor, photo_ion_block_references,
                  photo_ion_index, t_rad, w, phi_ik, t_electrons,
                  boltzmann_factor_photo_ion):
        # Used for initialization
        if alpha_stim_estimator is None:
            alpha_stim = self.calculate_from_dilute_bb(
                photo_ion_cross_sections, photo_ion_block_references,
                photo_ion_index, t_rad, w, t_electrons,
                boltzmann_factor_photo_ion
            )
            alpha_stim *= phi_ik.loc[alpha_stim.index]
        else:
            alpha_stim = alpha_stim_estimator * photo_ion_norm_factor
        return alpha_stim

    @staticmethod
    def calculate_from_dilute_bb(photo_ion_cross_sections,
                                 photo_ion_block_references,
                                 photo_ion_index, t_rad, w, t_electrons,
                                 boltzmann_factor_photo_ion):
        nu = photo_ion_cross_sections['nu']
        x_sect = photo_ion_cross_sections['x_sect']
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        j_nus *= boltzmann_factor_photo_ion
        alpha_stim = j_nus.multiply(
            4. * np.pi * x_sect / nu / const.h.cgs.value, axis=0
        )
        alpha_stim = integrate_array_by_blocks(alpha_stim.values, nu.values,
                                               photo_ion_block_references)
        alpha_stim = pd.DataFrame(alpha_stim, index=photo_ion_index)
        return alpha_stim


class BaseRecombTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_recomb : Pandas DataFrame, dtype float
               The unnormalized transition probabilities for
               spontaneous recombination.
    """
    outputs = ('p_recomb', )
    transition_probabilities_outputs = ('p_recomb', )
    latex_name = ('p^{\\textrm{recomb}}', '')

    def calculate(self, alpha_sp, nu_i, energy_i, photo_ion_idx):
        p_recomb_deac = alpha_sp.multiply(nu_i, axis=0) * const.h.cgs.value
        p_recomb_deac = self.set_index(p_recomb_deac, photo_ion_idx,
                                       transition_type=-1)

        p_recomb_internal = alpha_sp.multiply(energy_i, axis=0)
        p_recomb_internal = self.set_index(p_recomb_internal, photo_ion_idx,
                                           transition_type=0)
        p_recomb = pd.concat([p_recomb_deac, p_recomb_internal])
        return p_recomb


class BasePhotoIonTransProbs(TransitionProbabilitiesProperty,
                             IndexSetterMixin):
    """
    Attributes
    ----------
    p_photo_ion : Pandas DataFrame, dtype float
                  The unnormalized transition probabilities for
                  radiative ionization.
    """
    outputs = ('p_photo_ion', )
    transition_probabilities_outputs = ('p_photo_ion', )
    latex_name = ('p^{\\textrm{photo_ion}}', )

    def calculate(self, gamma_corr, nu_i, photo_ion_idx):
        p_photo_ion = gamma_corr.multiply(nu_i, axis=0) * const.h.cgs.value
        p_photo_ion = self.set_index(p_photo_ion, photo_ion_idx, reverse=False)
        return p_photo_ion


class CorrPhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma_corr : Pandas DataFrame, dtype float
        The rate coefficient for radiative ionization corrected for
        stimulated recombination.
    """
    outputs = ('gamma_corr',)
    latex_name = ('\\gamma_\\mathrm{corr}',)

    def calculate(self, gamma, alpha_stim, electron_densities,
                  ion_number_density, level_number_density):
        n_k_index = get_ion_multi_index(alpha_stim.index)
        n_k = ion_number_density.loc[n_k_index].values
        n_i = level_number_density.loc[alpha_stim.index].values
        gamma_corr = gamma - alpha_stim * n_k * electron_densities / n_i
        num_neg_elements = (gamma_corr < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException(
                "Negative values in CorrPhotoIonRateCoeff."
            )
        return gamma_corr


class PhotoIonEstimatorsNormFactor(ProcessingPlasmaProperty):
    outputs = ('photo_ion_norm_factor',)
    latex = ('\\frac{1}}}{'
             'time_\\textrm{simulation} volume h}')

    @staticmethod
    def calculate(time_simulation, volume):
        return (time_simulation * volume * const.h.cgs.value)**-1


class PhotoIonRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    gamma_estimator : Pandas DataFrame, dtype float
                      Unnormalized MC estimator for the rate coefficient
                      for radiative ionization.
    """
    outputs = ('gamma_estimator',)


class StimRecombRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    alpha_stim_estimator : Pandas DataFrame, dtype float
                           Unnormalized MC estimator for the rate coefficient
                           for stimulated recombination.
    """
    outputs = ('alpha_stim_estimator',)


class BfHeatingRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    bf_heating_coeff_estimator : Pandas DataFrame, dtype float
                                 Unnormalized MC estimator for the rate
                                 coefficient for bound-free heating.
    """
    outputs = ('bf_heating_coeff_estimator',)


class CollExcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_exc_coeff : Pandas DataFrame, dtype float
        Rate coefficient for collisional excitation.
    """
    outputs = ('coll_exc_coeff',)
    latex_name = ('c_{lu}',)

    def calculate(self, yg_interp, yg_index, t_electrons, delta_E_yg):
        yg = yg_interp(t_electrons)
        k_B = const.k_B.cgs.value
        boltzmann_factor = np.exp(
            - delta_E_yg.values[np.newaxis].T / (t_electrons * k_B)
        )
        q_ij = 8.629e-6 / np.sqrt(t_electrons) * yg * boltzmann_factor  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181
        return pd.DataFrame(q_ij, index=yg_index)


class CollDeexcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_deexc_coeff : Pandas DataFrame, dtype float
        Rate coefficient for collisional deexcitation.
    """
    outputs = ('coll_deexc_coeff',)
    latex_name = ('c_{ul}',)

    def calculate(self, thermal_lte_level_boltzmann_factor, coll_exc_coeff):
        ll_index = coll_exc_coeff.index.droplevel('level_number_upper')
        lu_index = coll_exc_coeff.index.droplevel('level_number_lower')

        n_lower_prop = thermal_lte_level_boltzmann_factor.loc[ll_index].values
        n_upper_prop = thermal_lte_level_boltzmann_factor.loc[lu_index].values

        coll_deexc_coeff = coll_exc_coeff * n_lower_prop / n_upper_prop
        return coll_deexc_coeff


class BaseCollisionTransProbs(TransitionProbabilitiesProperty,
                              IndexSetterMixin):
    """
    Attributes
    ----------
    p_coll : Pandas DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional excitation.
    """
    outputs = ('p_coll', )
    transition_probabilities_outputs = ('p_coll', )
    latex_name = ('p^{\\textrm{coll}}', )

    def calculate(self, coll_exc_coeff, coll_deexc_coeff, yg_idx,
                  electron_densities, delta_E_yg, atomic_data,
                  level_number_density):
        p_deexc_deac = (coll_deexc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0)
        p_deexc_deac = self.set_index(p_deexc_deac, yg_idx)
        p_deexc_deac = p_deexc_deac.groupby(level=[0]).sum()
        index_dd = pd.MultiIndex.from_product(
            [p_deexc_deac.index.values, ['k'], [0]],
            names=list(yg_idx.columns) + ['transition_type']
        )
        p_deexc_deac = p_deexc_deac.set_index(index_dd)

        ll_index = coll_deexc_coeff.index.droplevel('level_number_upper')
        energy_lower = atomic_data.levels.energy.loc[ll_index]
        p_deexc_internal = (coll_deexc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0)
        p_deexc_internal = self.set_index(p_deexc_internal, yg_idx,
                                          transition_type=0, reverse=True)

        p_exc_internal = (coll_exc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0)
        p_exc_internal = self.set_index(p_exc_internal, yg_idx,
                                        transition_type=0, reverse=False)
        p_exc_cool = (coll_exc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0)
        p_exc_cool = p_exc_cool * level_number_density.loc[ll_index].values
        p_exc_cool = self.set_index(p_exc_cool, yg_idx, reverse=False)
        p_exc_cool = p_exc_cool.groupby(level='destination_level_idx').sum()
        exc_cool_index = pd.MultiIndex.from_product(
            [['k'], p_exc_cool.index.values, [0]],
            names=list(yg_idx.columns) + ['transition_type']
        )
        p_exc_cool = p_exc_cool.set_index(exc_cool_index)
        p_coll = pd.concat(
            [p_deexc_deac, p_deexc_internal, p_exc_internal, p_exc_cool]
        )
        return p_coll


class AdiabaticCoolingRate(TransitionProbabilitiesProperty):
    """
    Attributes
    ----------
    C_adiabatic : Pandas DataFrame, dtype float
        The adiabatic cooling rate of the electron gas.
    """
    outputs = ('C_adiabatic', )
    transition_probabilities_outputs = ('C_adiabatic', )
    latex_name = ('C_{\\textrm{adiabatic}}', )

    def calculate(self, electron_densities, t_electrons, time_explosion):
        C_adiabatic = (3. * electron_densities * const.k_B.cgs.value *
                       t_electrons) / time_explosion

        C_adiabatic = cooling_rate_series2dataframe(
            C_adiabatic, destination_level_idx='adiabatic'
        )
        return C_adiabatic


class FreeFreeCoolingRate(TransitionProbabilitiesProperty):
    """
    Attributes
    ----------
    C_ff : Pandas DataFrame, dtype float
        The free-free cooling rate of the electron gas.
    """
    outputs = ('C_ff', )
    transition_probabilities_outputs = ('C_ff', )
    latex_name = ('C^{\\textrm{ff}}', )

    def calculate(self, ion_number_density, electron_densities,
                  t_electrons):
        ff_cooling_factor = self._calculate_ff_cooling_factor(
            ion_number_density, electron_densities
        )
        C_ff = 1.426e-27 * np.sqrt(t_electrons) * ff_cooling_factor
        C_ff = cooling_rate_series2dataframe(
            C_ff, destination_level_idx='ff'
        )
        return C_ff

    @staticmethod
    def _calculate_ff_cooling_factor(ion_number_density, electron_densities):
        ion_charge = ion_number_density.index.get_level_values(1).values
        factor = electron_densities * ion_number_density.multiply(
            ion_charge ** 2, axis=0).sum()
        return factor


class BoundFreeOpacity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    chi_bf : Pandas DataFrame, dtype float
    """
    outputs = ('chi_bf',)
    latex_name = ('\\chi^{\\textrm{bf}}',)

    def calculate(self, photo_ion_cross_sections, t_electrons,
                  phi_ik, level_number_density, lte_level_number_density,
                  boltzmann_factor_photo_ion):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values

        n_i = level_number_density.loc[photo_ion_cross_sections.index]
        lte_n_i = lte_level_number_density.loc[photo_ion_cross_sections.index]
        chi_bf = (n_i - lte_n_i * boltzmann_factor_photo_ion).multiply(
            x_sect, axis=0
        )

        num_neg_elements = (chi_bf < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException(
                "Negative values in bound-free opacity."
            )
        return chi_bf


class LevelNumberDensityLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    lte_level_number_density : Pandas DataFrame, dtype float
    """
    outputs = ('lte_level_number_density',)
    latex_name = ('n_{\\textrm{i}}^*',)

    # TODO: only do this for continuum species
    def calculate(self, electron_densities, phi_ik, ion_number_density):
        next_higher_ion_index = get_ion_multi_index(
            phi_ik.index, next_higher=True
        )
        # TODO: Check that n_k is correct (and not n_k*)
        lte_level_number_density = (
            phi_ik * ion_number_density.loc[next_higher_ion_index].values
        ).multiply(electron_densities, axis=1)
        return lte_level_number_density


class PhotoIonBoltzmannFactor(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    boltzmann_factor_photo_ion : Pandas DataFrame, dtype float
    """
    outputs = ('boltzmann_factor_photo_ion',)

    def calculate(self, photo_ion_cross_sections, t_electrons):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values

        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        return boltzmann_factor
