import numpy as np
import pytest

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_level_number_density_lte(level_number_density, ion_number_density):
    assert np.allclose(
        level_number_density.loc[2].loc[0].loc[0],
        ion_number_density.loc[2].loc[0],
    )
    assert np.allclose(level_number_density.sum(), ion_number_density.sum())
