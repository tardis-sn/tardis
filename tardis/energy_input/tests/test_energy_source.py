import pytest
import numpy as np
import numpy.testing as npt
import pandas as pd
import os

from tardis.energy_input.energy_source import (
    read_nuclear_dataframe,
    get_type_property,
    create_energy_cdf,
    sample_energy_distribution,
    setup_input_energy,
)


def test_read_nuclear_dataframe(tardis_ref_path):
    """
    Parameters
    ----------
    path : str
    """
    path = os.path.join(tardis_ref_path, "nuclear_data/simple_nuclear.h5")
    actual = read_nuclear_dataframe(path)
    expected = pd.read_hdf(path, key="decay_radiation")
    pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    ["type_of_radiation", "property", "expected"],
    [
        (
            "'gamma_rays'",
            "energy",
            3177.28,
        ),
        (
            "'e+'",
            "energy",
            65.0,
        ),
        (
            "'gamma_rays'",
            "intensity",
            0.0111,
        ),
    ],
)
def test_get_type_property(
    tardis_ref_path, type_of_radiation, property, expected
):
    """
    Parameters
    ----------
    path : str
    type_of_radiation : str
    property : str
    expected : float
    """
    path = os.path.join(tardis_ref_path, "nuclear_data/simple_nuclear.h5")
    nuclear_df = pd.read_hdf(path, key="decay_radiation")
    actual = get_type_property(nuclear_df, type_of_radiation, property)

    npt.assert_equal(actual[0], expected)


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
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_setup_input_energy():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_intensity_ratio():
    assert False
