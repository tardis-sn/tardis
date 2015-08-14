import logging

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelNumberDensity', 'LevelNumberDensityHeNLTE']

class LevelNumberDensity(ProcessingPlasmaProperty):
    """
    Outputs:
    level_number_density : Pandas DataFrame
    """
    outputs = ('level_number_density',)
    latex_name = ('N_{i,j,k}',)
    latex_formula = ('N_{i,j}\\dfrac{bf_{i,j,k}}{Z_{i,j}}',)

    def calculate(self, level_boltzmann_factor, ion_number_density,
        levels, partition_function):
        partition_function_broadcast = partition_function.ix[
            levels.droplevel(2)].values
        level_population_fraction = level_boltzmann_factor /\
            partition_function_broadcast
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast

class LevelNumberDensityHeNLTE(ProcessingPlasmaProperty):
    """
    Outputs:
    level_number_density : Pandas DataFrame
    """
    outputs = ('level_number_density',)
    latex_name = ('N_{i,j,k}',)
    latex_formula = ('N_{i,j}\\dfrac{bf_{i,j,k}}{Z_{i,j}}',)

    def calculate(self, level_boltzmann_factor, ion_number_density,
        levels, partition_function, helium_population):
        partition_function_broadcast = partition_function.ix[
            levels.droplevel(2)].values
        level_population_fraction = level_boltzmann_factor /\
            partition_function_broadcast
        level_population_fraction.ix[2].update(helium_population)
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast