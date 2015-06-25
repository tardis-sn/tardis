import numpy as np

def test_level_population(level_population):
    assert np.allclose(level_population.ix[2].ix[0].ix[0], 1.0)
    assert np.allclose(level_population.ix[2].ix[0].sum(), 1.0)

def test_level_number_density_lte(level_number_density,
    ion_number_density):
    assert np.allclose(level_number_density.ix[2].ix[0].ix[0],
        ion_number_density.ix[2].ix[0])
    assert np.allclose(level_number_density.sum(),
        ion_number_density.sum())
