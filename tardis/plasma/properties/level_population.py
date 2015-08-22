import logging

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

    def calculate():
        pass

    def __init__(self, plasma_parent, helium_treatment='dilute-lte'):
        super(LevelNumberDensity, self).__init__(plasma_parent)
        if hasattr(self.plasma_parent, 'plasma_properties_dict'):
            if 'HeliumNLTE' in \
                self.plasma_parent.plasma_properties_dict.keys():
                    helium_treatment=='recomb-nlte'
        if helium_treatment=='recomb-nlte' or 'numerical-nlte':
            self.calculate = self._calculate_helium_nlte
        elif helium_treatment=='dilute-lte':
            self.calculate = self._calculate_dilute_lte
        self._update_inputs()

    def _calculate_dilute_lte(self, level_boltzmann_factor, ion_number_density,
        levels, partition_function):
        partition_function_broadcast = partition_function.ix[
            levels.droplevel(2)].values
        level_population_fraction = level_boltzmann_factor /\
            partition_function_broadcast
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast

    def _calculate_helium_nlte(self, level_boltzmann_factor,
        ion_number_density, levels, partition_function, helium_population):
        level_number_density = self._calculate_dilute_lte(
            level_boltzmann_factor, ion_number_density, levels,
            partition_function)
        if helium_population is not None:
            level_number_density.ix[2].update(helium_population)
        return level_number_density