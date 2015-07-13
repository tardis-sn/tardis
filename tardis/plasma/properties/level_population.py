import logging

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelPopulationFraction', 'LevelNumberDensity']

class LevelPopulationFraction(ProcessingPlasmaProperty):
    """
    Outputs:
    level_population_fraction : Pandas DataFrame
    """
    outputs = ('level_population_fraction',)
    latexformula = (r'$N_{i, j, k} = \frac{g_{i, j, k} '
                         r'e^{-\beta_\textrm{rad} E_{i, j, k}}}{Z_{i, j}}$')

    def calculate(self, levels, partition_function, level_boltzmann_factor):
        partition_function_broadcast = partition_function.ix[
            levels.index.droplevel(2)].values
        return level_boltzmann_factor / partition_function_broadcast

class LevelNumberDensity(ProcessingPlasmaProperty):
    """
    Outputs:
    level_number_density : Pandas DataFrame
    """
    outputs = ('level_number_density',)

    def calculate(self, level_population_fraction, ion_number_density):
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast