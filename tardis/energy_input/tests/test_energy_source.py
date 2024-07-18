import numpy as np
import numpy.testing as npt
import pytest

from tardis.energy_input.samplers import (
    create_energy_cdf,
)


@pytest.mark.parametrize(
    ["energy", "intensity", "expected_cdf"],
    [
        (np.array([100.0, 50.0]), np.array([1.0, 1.0]), np.array([0.5, 1.0])),
        (np.array([50.0, 100.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0])),
    ],
)
def test_create_energy_cdf(energy, intensity, expected_cdf):
    """
    Parameters
    ----------
    energy : One-dimensional Numpy Array, dtype float
    intensity : One-dimensional Numpy Array, dtype float
    expected_cdf : One-dimensional Numpy Array, dtype float
    """
    actual_energy, actual_cdf = create_energy_cdf(energy, intensity)
    expected_energy = np.sort(energy)

    npt.assert_array_almost_equal_nulp(actual_cdf, expected_cdf)
    npt.assert_array_almost_equal_nulp(actual_energy, expected_energy)


@pytest.mark.xfail(reason="To be implemented")
def test_sample_energy_distribution():
    """To test the energy distribution sample"""
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_setup_input_energy():
    """To test setting up the input energy"""
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_intensity_ratio():
    """To test the intensity ratio"""
    assert False
