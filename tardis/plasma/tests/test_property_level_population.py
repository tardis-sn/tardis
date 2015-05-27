import numpy as np

from tardis.plasma.properties.level_population import (LevelPopulationLTE,
LevelNumberDensity)

def test_level_population_lte(level_population_lte):
    assert np.allclose(level_population_lte.ix[2].ix[0].ix[0], 1.0)
    assert np.allclose(level_population_lte.ix[2].ix[0].sum(), 1.0)

def test_level_number_density(level_population_lte, ion_number_density):
    level_number_density_module = LevelNumberDensity(None)
    level_number_density = level_number_density_module.calculate(
        level_population_lte, ion_number_density)
    assert np.allclose(level_number_density.ix[2].ix[0].ix[0],
        ion_number_density.ix[2].ix[0])
    assert np.allclose(level_number_density.sum(), ion_number_density.sum())
