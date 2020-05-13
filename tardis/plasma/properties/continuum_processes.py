import logging

import numpy as np
import pandas as pd

from numba import prange, njit
from astropy import constants as const

from tardis.plasma.exceptions import PlasmaException
from tardis.plasma.properties.base import (ProcessingPlasmaProperty,
                                           Input)
from tardis.plasma.properties.j_blues import JBluesDiluteBlackBody

__all__ = ['SpontRecombRateCoeff', 'StimRecombRateCoeff', 'PhotoIonRateCoeff',
           'PhotoIonEstimatorsNormFactor', 'PhotoIonRateCoeffEstimator',
           'StimRecombRateCoeffEstimator', 'CorrPhotoIonRateCoeff']

logger = logging.getLogger(__name__)

njit_dict = {'fastmath': False, 'parallel': False}


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrates a function f defined at locations x over blocks
    given in block_references.

    Parameters
    ----------
    f : Two-dimensional Numpy Array, dtype float
    x : One-dimensional Numpy Array, dtype float
    block_references : One-dimensional Numpy Array, dtype int

    Returns
    -------
    integrated : Two-dimensional Numpy Array, dtype float

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
    Integrates a function f defined at locations x over blocks
    given in block_references.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)
    next_higher : bool
        If true use ion number of next higher ion, else use ion_number from
        multi_index_full.

    Returns
    -------
    multi_index : Pandas MultiIndex (atomic_number, ion_number)

    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


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
                  photo_ion_block_references, photo_ion_index, phi_ik):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values

        alpha_sp = (8 * np.pi * x_sect * nu ** 2 / (const.c.cgs.value) ** 2)
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        alpha_sp = alpha_sp * boltzmann_factor
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
                  photo_ion_index, t_rad, w, phi_ik, t_electrons):
        # Used for initialization
        if alpha_stim_estimator is None:
            alpha_stim = self.calculate_from_dilute_bb(
                photo_ion_cross_sections, photo_ion_block_references,
                photo_ion_index, t_rad, w, t_electrons
            )
            alpha_stim *= phi_ik.loc[alpha_stim.index]
        else:
            alpha_stim = alpha_stim_estimator * photo_ion_norm_factor
        return alpha_stim

    @staticmethod
    def calculate_from_dilute_bb(photo_ion_cross_sections,
                                 photo_ion_block_references,
                                 photo_ion_index, t_rad, w, t_electrons):
        nu = photo_ion_cross_sections['nu']
        x_sect = photo_ion_cross_sections['x_sect']
        boltzmann_factor = np.exp(-nu.values[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        j_nus *= boltzmann_factor
        alpha_stim = j_nus.multiply(
            4. * np.pi * x_sect / nu / const.h.cgs.value, axis=0
        )
        alpha_stim = integrate_array_by_blocks(alpha_stim.values, nu.values,
                                               photo_ion_block_references)
        alpha_stim = pd.DataFrame(alpha_stim, index=photo_ion_index)
        return alpha_stim


class CorrPhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
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
