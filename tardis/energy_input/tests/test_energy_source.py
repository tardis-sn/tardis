import pytest
import numpy as np
import numpy.testing as npt

from tardis.energy_input.energy_source import (
    read_nuclear_dataframe,
    get_type_property,
    create_energy_cdf,
    sample_energy_distribution,
    setup_gamma_ray_energy,
)


@pytest.mark.xfail(reason="To be implemented")
def test_read_nuclear_dataframe():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_get_type_property():
    assert False


@pytest.mark.parametrize(
    ["energy", "intensity", "expected_cdf"],
    [
        (np.array([100.0, 50.0]), np.array([1.0, 1.0]), np.array([0.5, 1.0])),
        (np.array([50.0, 100.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0])),
    ],
)
def test_create_energy_cdf(energy, intensity, expected_cdf):
    actual_energy, actual_cdf = create_energy_cdf(energy, intensity)
    expected_energy = np.sort(energy)

    npt.assert_array_almost_equal_nulp(actual_cdf, expected_cdf)
    npt.assert_array_almost_equal_nulp(actual_energy, expected_energy)


@pytest.mark.xfail(reason="To be implemented")
def test_sample_energy_distribution():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_setup_gamma_ray_energy():
    assert False