from abc import ABCMeta, abstractmethod
import logging

from astropy import constants as const
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

class BasePlasmaProperty(object):
    __metaclass__ = ABCMeta

    def __init__(self, plasma_parent):
        self.plasma_parent = plasma_parent

    def update(self):
        args = [getattr(self.plasma_parent, item) for item in self.inputs]
        self.value = self.calculate(*args)

    @abstractmethod
    def calculate(self, *args, **kwargs):
        raise NotImplementedError('This method needs to be implemented by ')

    def get_label(self):
        return "Name: {0}\nType: {1}\n{2}".format(self.name.replace('_', r'\\_'),
                                                  self.type_str,
                                                  getattr(self,
                                                          'latex_str', ''))

class BetaRadiation(BasePlasmaProperty):

    name = 'beta_rad'
    inputs = ['t_rad']
    type_str = 'numpy.array'
    latex_str = '$\\frac{1}{K_B T_\\textrm{rad}}$'

    def __init__(self, plasma_parent):
        super(BetaRadiation, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value


    def calculate(self, t_rad):
        return (1 / (self.k_B_cgs * t_rad))

class GElectron(BasePlasmaProperty):

    name = 'g_electron'
    inputs = ['beta_rad']

    @staticmethod
    def calculate(self, beta_rad):
        return ((2 * np.pi * const.m_e.cgs.value / beta_rad) /
                (const.h.cgs.value ** 2)) ** 1.5

class LevelBoltzmannFactor(BasePlasmaProperty):
    """
    Calculate the level population Boltzmann factor
    """

    name = 'level_boltzmann_factor'
    inputs = ['levels', 'beta_rad']
    type_str = 'pandas.DataFrame'

    def calculate(self, levels, beta_rad):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        return level_boltzmann_factor

class LTEPartitionFunction(BasePlasmaProperty):
    name = 'function'
    inputs = ['levels', 'level_boltzmann_factor']
    type_str = 'NA'

    @staticmethod
    def calculate(levels, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=['atomic_number', 'ion_number']).sum()

class DiluteLTEPartitionFunction(BasePlasmaProperty):
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
    inputs = ['levels', 'level_boltzmann_factor']
    type_str = 'NA'

    @staticmethod
    def calculate(w, levels, level_boltzmann_factor):
        metastable = levels.metastable
        partition_functions = level_boltzmann_factor[metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()

        partition_functions_non_meta = w * level_boltzmann_factor[~metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()

        partition_functions.ix[
            partition_functions_non_meta.index] += partition_functions_non_meta

        #if self.nlte_config is not None   and self.nlte_config.species != [] and not initialize_nlte:
        #    for species in self.nlte_config.species:
        #        partition_functions.ix[species] = self.atom_data.levels.g.ix[species].ix[0] * \
        #                                               (self.level_populations.ix[species] /
        #                                                self.level_populations.ix[species].ix[0]).sum()


class PartitionFunction(BasePlasmaProperty):
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



class NumberDensity(BasePlasmaProperty):
    name = 'number_density'
    inputs = ['atomic_mass', 'abundance']
    type_str = 'pandas.DataFrame'

    def calculate(self, atomic_mass, abundance):
        pass

class SelectedAtoms(BasePlasmaProperty):
    name = 'selected_atoms'
    inputs = ['abundance']
    type_str = 'iterable'

    def calculate(self, abundance):
        return self.plasma_parent.abundance.index