import os
import pandas as pd
import pytest
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from tardis.io.config_reader import Configuration
from tardis.model.density import HomologousDensity


###
# Save and Load
###

def data_path(filename):
    return os.path.join('tardis/io/tests/data/', filename)

@pytest.fixture(scope="module")
def hdf_file_path(tmpdir_factory):
    path = tmpdir_factory.mktemp('hdf_buffer').join('density.hdf')
    return str(path)

@pytest.fixture(scope="module")
def homologous_density():
    filename = 'tardis_configv1_verysimple.yml'
    config = Configuration.from_yaml(data_path(filename)) 
    density =  HomologousDensity.from_config(config)
    return density

@pytest.fixture(scope="module",autouse=True)
def to_hdf_buffer(hdf_file_path,homologous_density):
    homologous_density.to_hdf(hdf_file_path)

def test_hdf_density_0(hdf_file_path, homologous_density):
    actual = homologous_density.density_0
    if hasattr(actual, 'cgs'):
        actual = actual.cgs.value
    path = os.path.join('homologous_density','density_0')
    expected = pd.read_hdf(hdf_file_path, path)
    assert_array_almost_equal(actual, expected.values)

def test_hdf_time_0(hdf_file_path, homologous_density):
    actual = homologous_density.time_0
    if hasattr(actual, 'cgs'):
        actual = actual.cgs.value
    path = os.path.join('homologous_density','scalars')
    expected = pd.read_hdf(hdf_file_path, path)['time_0']
    assert_almost_equal(actual, expected)