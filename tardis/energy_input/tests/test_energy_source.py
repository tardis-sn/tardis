import pytest
import numpy as np
import numpy.testing as npt
import pandas as pd

from tardis.energy_input.energy_source import (
    read_nuclear_dataframe,
    get_type_property,
    create_energy_cdf,
    sample_energy_distribution,
    setup_gamma_ray_energy,
)

# TODO make dummy nuclear_data


@pytest.mark.xfail(reason="To be implemented")
def test_read_nuclear_dataframe(path):
    actual = read_nuclear_dataframe(path)
    expected = pd.read_hdf(path, key="decay_radiation")


@pytest.mark.parametrize(
    ["path", "type_of_radiation", "property"],
    [
        (
            "/home/afullard/Downloads/tardisnuclear/decay_radiation.h5",
            "'gamma_rays'",
            "energy",
        ),
        (
            "/home/afullard/Downloads/tardisnuclear/decay_radiation.h5",
            "'e+'",
            "energy",
        ),
        (
            "/home/afullard/Downloads/tardisnuclear/decay_radiation.h5",
            "'gamma_rays'",
            "intensity",
        ),
    ],
)
def test_get_type_property(path, type_of_radiation, property):
    nuclear_df = pd.read_hdf(path, key="decay_radiation")

    expected = nuclear_df.query("type==" + type_of_radiation)[property].values
    actual = get_type_property(nuclear_df, type_of_radiation, property)

    npt.assert_array_equal(actual, expected)


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