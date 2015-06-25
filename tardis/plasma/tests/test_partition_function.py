import numpy as np

def test_level_boltzmann_factor_lte(level_boltzmann_factor_lte, levels):
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[0],
        levels.ix[2].ix[0].ix[0][1])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[1].ix[0],
        levels.ix[2].ix[1].ix[0][1])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[10],
        7.6218092841240705e-12)

def test_level_boltzmann_factor_dilute_lte(level_boltzmann_factor_dilute_lte,
    levels):
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[0],
        levels.ix[2].ix[0].ix[0][1])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[1].ix[0],
        levels.ix[2].ix[1].ix[0][1])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[10],
        3.8109046420620353e-12)

def test_lte_partition_function(partition_function, levels):
    assert np.allclose(partition_function.ix[2].ix[0], 1.0)
    assert np.allclose(partition_function.ix[2].ix[1], 2.0)
    assert np.allclose(partition_function.ix[2].ix[2], 1.0)
