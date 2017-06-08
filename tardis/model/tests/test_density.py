import os

import pytest
from astropy.tests.helper import assert_quantity_allclose

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
    homologous_density.to_hdf(hdf_file_path,'')


@pytest.fixture(scope="module")
def from_hdf_buffer(hdf_file_path):
    hdf_buffer = HomologousDensity.from_hdf('',hdf_file_path)
    return hdf_buffer


homologous_density_attrs = ['density_0', 'time_0']

@pytest.mark.parametrize("attr", homologous_density_attrs)
def test_hdf_homologous_density(from_hdf_buffer, homologous_density, attr):
    if hasattr(homologous_density, attr):
        assert_quantity_allclose(getattr(homologous_density, attr), getattr(
            from_hdf_buffer, attr))
