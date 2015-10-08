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

logger = logging.getLogger(__name__)


class PhotoIonizationRates(object):
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
        # corrected_photoion_coeff = corrected_photoion_coeff.reset_index()\
        #    .groupby(('atomic_number','ion_number','level_number'))
        corrected_photoion_coeff = corrected_photoion_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(self.t_rads)):
            tmp[i] = corrected_photoion_coeff.apply(lambda sub: simps(sub[i], sub['nu']))
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