import numpy as np

from tardis.plasma.properties import LevelBoltzmannFactor, LTEPartitionFunction


#def test_level_boltzmann_factor1(standard_lte_plasma_he_db):
#    level_boltzmann_factor_module = LevelBoltzmannFactor(None)
#    levels = standard_lte_plasma_he_db.levels
#    factor = level_boltzmann_factor_module.calculate(levels, standard_lte_plasma_he_db.beta_rad)
#    assert np.allclose(factor.loc[2, 0, 0], levels.loc[2, 0, 0].g)
#    assert np.allclose(factor.loc[2, 1, 0], levels.loc[2, 1, 0].g)
#    assert np.allclose(factor.loc[2, 0, 10], 7.6218092841240705e-12)

#def test_lte_partition_function1(standard_lte_plasma_he_db):
#    partition_function = LTEPartitionFunction.calculate(
#        standard_lte_plasma_he_db.levels,
#        standard_lte_plasma_he_db.level_boltzmann_factor)
#    assert 0
