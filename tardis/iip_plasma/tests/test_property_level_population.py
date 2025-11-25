import numpy as np


def test_level_number_density_lte(level_number_density, ion_number_density):
    assert np.allclose(
        level_number_density.iloc[2].iloc[0].iloc[0],
        ion_number_density.iloc[2].iloc[0],
    )
    assert np.allclose(level_number_density.sum(), ion_number_density.sum())
