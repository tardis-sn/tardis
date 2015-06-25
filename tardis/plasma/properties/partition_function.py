import logging

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelBoltzmannFactorLTE', 'LevelBoltzmannFactorDiluteLTE',
           'PartitionFunction']

class LevelBoltzmannFactorLTE(ProcessingPlasmaProperty):
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

class LevelBoltzmannFactorDiluteLTE(ProcessingPlasmaProperty):

    name = 'level_boltzmann_factor'

    @staticmethod
    def calculate(levels, beta_rad, w):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        level_boltzmann_factor[~levels.metastable] *= w
        return level_boltzmann_factor

class PartitionFunction(ProcessingPlasmaProperty):
    name = 'partition_function'
    latex_name = '$Z_{i, j}$'

    latex_formula = (r'$Z_{i, j} = \sum_{k=1}^n g_{i, j, k} '
                     r'e^{E_{i, j, k} \times \beta_\textrm{rad}}$')

    @staticmethod
    def calculate(levels, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=['atomic_number', 'ion_number']).sum()
