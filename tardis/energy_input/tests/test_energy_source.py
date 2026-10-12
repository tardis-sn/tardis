import numpy as np
import numpy.testing as npt
import pytest

from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.energy_input.samplers import (
    PositroniumSampler,
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
    expected_energy = np.sort(energy, kind=SORTING_ALGORITHM)

    npt.assert_array_almost_equal_nulp(actual_cdf, expected_cdf)
    npt.assert_array_almost_equal_nulp(actual_energy, expected_energy)


def test_positronium_sampler_energy_weights_photon_pdf() -> None:
    sampler = PositroniumSampler(n_grid=1000)
    energy_pdf = sampler.x_grid * sampler.pdf(sampler.x_grid)
    expected_cdf = np.cumsum(energy_pdf)
    expected_cdf /= expected_cdf[-1]

    npt.assert_allclose(sampler.cdf_grid, expected_cdf)
