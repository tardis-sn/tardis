import os
import pandas as pd
import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.montecarlo import MontecarloRunner


###
# Save and Load
###

def data_path(filename):
    return os.path.abspath(os.path.join('tardis/io/tests/data/', filename))


@pytest.fixture(scope="module")
def hdf_file_path(tmpdir_factory):
    path = tmpdir_factory.mktemp('hdf_buffer').join('runner.hdf')
    return str(path)


@pytest.fixture(scope="module")
def config():
    filename = 'tardis_configv1_verysimple.yml'
    return Configuration.from_yaml(data_path(filename))


@pytest.fixture(scope="module")
def model(config):
    return Radial1DModel.from_config(config)


class Mock_plasma():
        def __init__(self):
            self.tau_sobolevs = np.zeros((29224, 20))


@pytest.fixture(scope="module")
def actual(config, model):
    runner = MontecarloRunner.from_config(config)
    runner.time_of_simulation = runner.calculate_time_of_simulation(model)
    runner.volume = model.volume
    plasma = Mock_plasma()
    runner._initialize_estimator_arrays(runner.volume.shape[0],
                                        plasma.tau_sobolevs.shape)
    runner._initialize_geometry_arrays(model)
    runner._initialize_packets(model.t_inner.value, no_of_packets=20)
    return runner

@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path, actual):
    actual.to_hdf(hdf_file_path, name='runner')

runner_properties = ['output_nu', 'output_energy', 'nu_bar_estimator',
                     'j_estimator', 'montecarlo_virtual_luminosity',
                     'last_interaction_in_nu',
                     'last_line_interaction_in_id',
                     'last_line_interaction_out_id',
                     'last_line_interaction_shell_id',
                     'packet_luminosity']

@pytest.mark.parametrize("attr", runner_properties)
def test_hdf_runner(hdf_file_path, actual, attr):
    actual_property = getattr(actual, attr)
    if hasattr(actual_property, 'cgs'):
        actual_property = actual_property.cgs.value
    path = os.path.join('runner', attr)
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual_property, expected.values)
