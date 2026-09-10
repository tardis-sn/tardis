import numpy as np
import pytest

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_level_boltzmann_factor_lte(level_boltzmann_factor_lte, levels):
    assert np.allclose(level_boltzmann_factor_lte.loc[2].loc[0].loc[0], 1)
    assert np.allclose(level_boltzmann_factor_lte.loc[2].loc[1].loc[0], 2)
    assert np.allclose(
        level_boltzmann_factor_lte.loc[2].loc[0].loc[10],
        7.6218092841240705e-12,
    )


def test_level_boltzmann_factor_dilute_lte(
    level_boltzmann_factor_dilute_lte, levels
):
    assert np.allclose(
        level_boltzmann_factor_dilute_lte.loc[2].loc[0].loc[0], 1
    )
    assert np.allclose(
        level_boltzmann_factor_dilute_lte.loc[2].loc[1].loc[0], 2
    )
    assert np.allclose(
        level_boltzmann_factor_dilute_lte.loc[2].loc[0].loc[10],
        3.8109046420620353e-12,
    )


def test_lte_partition_function(partition_function, levels):
    assert np.allclose(partition_function.loc[2].loc[0], 1.0)
    assert np.allclose(partition_function.loc[2].loc[1], 2.0)
    assert np.allclose(partition_function.loc[2].loc[2], 1.0)
