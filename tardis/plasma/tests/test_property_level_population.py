import numpy as np

def test_level_population_lte(level_population_lte):
    assert np.allclose(level_population_lte.ix[2].ix[0].ix[0], 1.0)
    assert np.allclose(level_population_lte.ix[2].ix[0].sum(), 1.0)

def test_level_number_density_lte(level_number_density_lte,
    ion_number_density_lte):
    assert np.allclose(level_number_density_lte.ix[2].ix[0].ix[0],
        ion_number_density_lte.ix[2].ix[0])
    assert np.allclose(level_number_density_lte.sum(),
        ion_number_density_lte.sum())
