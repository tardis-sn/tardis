import logging

import numpy as np
import pandas as pd

from tardis.plasma.base_properties import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)



class LevelBoltzmannFactor(ProcessingPlasmaProperty):
    """
    Calculate the level population Boltzmann factor

    .. math:
        {latex_formula}

    """

    name = 'level_boltzmann_factor'
    latex_formula = r'$g_{i, j, k} e^{E_{i, j, k} \times \beta_\textrm{rad}}$'

    @staticmethod
    def calculate(levels, beta_rad):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        return level_boltzmann_factor

class LTEPartitionFunction(ProcessingPlasmaProperty):
    name = 'partition_function'
    latex_name = '$Z_{i, j}$'

    latex_formula = (r'$Z_{i, j} = \sum_{k=1}^n g_{i, j, k} '
                     r'e^{E_{i, j, k} \times \beta_\textrm{rad}}$')

    @staticmethod
    def calculate(levels, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=['atomic_number', 'ion_number']).sum()

class DiluteLTEPartitionFunction(ProcessingPlasmaProperty):
    """
    Calculate partition functions for the ions using the following formula, where
    :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

    .. math::
        Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



    if self.initialize is True set the first time the partition functions are initialized.
    This will set a self.partition_functions and initialize with LTE conditions.


    Returns
    -------

    partition_functions : `~astropy.table.Table`
        with fields atomic_number, ion_number, partition_function

    """
    name = 'partition_function'
    inputs = ['w', 'levels', 'level_boltzmann_factor']
    latex_formula = (r'$Z_{i, j} = '
                     r'\begin{cases}'
                     r'\sum_{k=1}^n g_{i, j, k} e^{E_{i, j, k} '
                     r'\beta_\textrm{rad}} &'
                     r' \text{if level metastable} \\'
                     r'W \sum_{k=1}^n g_{i, j, k} e^{E_{i, j, k} '
                     r'\beta_\textrm{rad}} &'
                     r' \text{if level not metastable}'
                     r'\end{cases}$')

    @staticmethod
    def calculate(w, levels, level_boltzmann_factor):
        metastable = levels.metastable
        partition_functions = level_boltzmann_factor[metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()

        partition_functions_non_meta = w * level_boltzmann_factor[~metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()

        partition_functions.ix[
            partition_functions_non_meta.index] += partition_functions_non_meta


class PartitionFunction(ProcessingPlasmaProperty):
    """
    Calculate partition functions for the ions using the following formula, where
    :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

    .. math::
        Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



    if self.initialize is True set the first time the partition functions are initialized.
    This will set a self.partition_functions and initialize with LTE conditions.


    Returns
    -------

    partition_functions : `~astropy.table.Table`
        with fields atomic_number, ion_number, partition_function

    """

