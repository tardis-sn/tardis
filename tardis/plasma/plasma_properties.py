from astropy import constants as const
import numpy as np
import pandas as pd

class PlasmaProperty(object):
    def __init__(self, plasma_parent):
        self.plasma_parent = plasma_parent
        self.current_cycle_id = None

    def __call__(self, *args):
        if self.current_cycle_id == self.plasma_parent.current_cycle_id:
            return getattr(self, self.name)
        else:
            self.current_cycle_id = self.plasma_parent.current_cycle_id
            setattr(self, self.name, self.calculate(*args))



class BetaRadiation(PlasmaProperty):

    name = 'beta_rad'
    inputs = ['t_rad']

    def __init__(self, plasma_parent):
        super(BetaRadiation, self).__init__(self, plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value


    def calculate(self, t_rad):
        return (1 / (self.k_B_cgs * t_rad))

class LevelBoltzmannFactor(PlasmaProperty):
    """
    Calculate the level population Boltzmann factor
    """


    name = 'level_boltzmann_factor'
    inputs = ['levels', 'beta_rad']

    def calculate(self, levels, beta_rad):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        return level_boltzmann_factor

class PartitionFunction(PlasmaProperty):
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

    inputs = ['levels']


    def calculate(self, levels, ):

        level_population_proportional_array = levels.g.values[np.newaxis].T *\
                                              np.exp(np.outer(levels.energy.values, -self.beta_rads))
        level_population_proportionalities = pd.DataFrame(level_population_proportional_array,
                                                               index=self.atom_data.levels.index,
                                                               columns=np.arange(len(self.t_rads)), dtype=np.float64)


        #level_props = self.level_population_proportionalities

        partition_functions = level_population_proportionalities[self.atom_data.levels.metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()
        partition_functions_non_meta = self.ws * level_population_proportionalities[~self.atom_data.levels.metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()
        partition_functions.ix[partition_functions_non_meta.index] += partition_functions_non_meta
        if self.nlte_config is not None and self.nlte_config.species != [] and not initialize_nlte:
            for species in self.nlte_config.species:
                partition_functions.ix[species] = self.atom_data.levels.g.ix[species].ix[0] * \
                                                       (self.level_populations.ix[species] /
                                                        self.level_populations.ix[species].ix[0]).sum()

        return level_population_proportionalities, partition_functions


