import logging

import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelPopulation', 'LevelNumberDensity']


class LevelPopulation(ProcessingPlasmaProperty):

        outputs = ('level_population_fraction',)
        latexformula = (r'$N_{i, j, k} = \frac{g_{i, j, k} '
                         r'e^{-\beta_\textrm{rad} E_{i, j, k}}}{Z_{i, j}}$')

        @staticmethod
        def calculate(levels, partition_function, level_boltzmann_factor):

            partition_function_broadcast = partition_function.ix[
                levels.index.droplevel(2)].values

            return level_boltzmann_factor / partition_function_broadcast

class LevelNumberDensity(ProcessingPlasmaProperty):
    outputs = ('level_number_density',)

    @staticmethod
    def calculate(level_population_fraction, ion_number_density):
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast