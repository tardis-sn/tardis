import os
import pandas as pd
import pytest
from numpy.testing import assert_almost_equal

###
# Save and Load
###


@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path, simulation_verysimple):
    simulation_verysimple.simulation_state.to_hdf(hdf_file_path, overwrite=True)


def test_hdf_density_0(hdf_file_path, simulation_verysimple):
    actual = simulation_verysimple.simulation_state.density
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    path = "simulation_state/density"
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual, expected.values)


def test_hdf_time_0(hdf_file_path, simulation_verysimple):
    actual = simulation_verysimple.simulation_state.time_explosion
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    path = "simulation_state/scalars"
    expected = pd.read_hdf(hdf_file_path, path)["time_explosion"]
    assert_almost_equal(actual, expected)
