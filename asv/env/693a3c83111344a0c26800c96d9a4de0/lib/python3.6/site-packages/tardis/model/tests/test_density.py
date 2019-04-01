import os
import pandas as pd
import pytest
from numpy.testing import assert_almost_equal

###
# Save and Load
###

@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path,simulation_verysimple):
    simulation_verysimple.model.homologous_density.to_hdf(hdf_file_path)

def test_hdf_density_0(hdf_file_path, simulation_verysimple):
    actual = simulation_verysimple.model.homologous_density.density_0
    if hasattr(actual, 'cgs'):
        actual = actual.cgs.value
    path = os.path.join('homologous_density','density_0')
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual, expected.values)

def test_hdf_time_0(hdf_file_path, simulation_verysimple):
    actual = simulation_verysimple.model.homologous_density.time_0
    if hasattr(actual, 'cgs'):
        actual = actual.cgs.value
    path = os.path.join('homologous_density','scalars')
    expected = pd.read_hdf(hdf_file_path, path)['time_0']
    assert_almost_equal(actual, expected)