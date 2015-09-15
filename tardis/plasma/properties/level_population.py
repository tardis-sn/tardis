import logging
import pandas as pd
import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelNumberDensity']

class LevelNumberDensity(ProcessingPlasmaProperty):
    """
    Outputs:
    level_number_density : Pandas DataFrame
    """
    outputs = ('level_number_density',)
    latex_name = ('N_{i,j,k}',)
    latex_formula = ('N_{i,j}\\dfrac{bf_{i,j,k}}{Z_{i,j}}',)

    def calculate(self):
        pass

    def __init__(self, plasma_parent, helium_treatment='dilute-lte'):
        super(LevelNumberDensity, self).__init__(plasma_parent)
        if hasattr(self.plasma_parent, 'plasma_properties_dict'):
            if 'HeliumNLTE' in \
                self.plasma_parent.plasma_properties_dict.keys():
                    helium_treatment='recomb-nlte'
        if helium_treatment=='recomb-nlte':
            self.calculate = self._calculate_helium_nlte
        elif helium_treatment=='dilute-lte':
            self.calculate = self._calculate_dilute_lte
        self._update_inputs()

        self.initialize_indices = True


    def _initialize_indices(self, levels, partition_function):
        indexer = pd.Series(np.arange(partition_function.shape[0]),
                            index=partition_function.index)
        self._ion2level_idx = indexer.ix[levels.droplevel(2)].values

    def _calculate_dilute_lte(self, level_boltzmann_factor, ion_number_density,
        levels, partition_function):

        if self.initialize_indices:
            self._initialize_indices(levels, partition_function)
            self.initialize_indices = False

        partition_function_broadcast = partition_function.values[
            self._ion2level_idx]
        level_population_fraction = (level_boltzmann_factor.values
                                     / partition_function_broadcast)
        ion_number_density_broadcast = ion_number_density.values[
            self._ion2level_idx]
        level_number_density = (level_population_fraction *
                                ion_number_density_broadcast)
        return pd.DataFrame(level_number_density,
                            index=level_boltzmann_factor.index)

    def _calculate_helium_nlte(self, level_boltzmann_factor,
        ion_number_density, levels, partition_function, helium_population):
        level_number_density = self._calculate_dilute_lte(
            level_boltzmann_factor, ion_number_density, levels,
            partition_function)
        if helium_population is not None:
            level_number_density.ix[2].update(helium_population)
        return level_number_density