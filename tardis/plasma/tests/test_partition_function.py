import numpy as np

def test_level_boltzmann_factor(level_boltzmann_factor, levels):
    assert np.allclose(level_boltzmann_factor.ix[2].ix[0].ix[0],
        levels.ix[2].ix[0].ix[0][1])
    assert np.allclose(level_boltzmann_factor.ix[2].ix[1].ix[0],
        levels.ix[2].ix[1].ix[0][1])
    assert np.allclose(level_boltzmann_factor.ix[2].ix[0].ix[10],
        7.6218092841240705e-12)

def test_lte_partition_function(partition_function_lte, levels):
    assert np.allclose(partition_function_lte.ix[2].ix[0], 1.0)
    assert np.allclose(partition_function_lte.ix[2].ix[1], 2.0)
    assert np.allclose(partition_function_lte.ix[2].ix[2], 1.0)
