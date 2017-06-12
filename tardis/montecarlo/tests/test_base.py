import os

import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

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
    actual.to_hdf(hdf_file_path, 'runner')


@pytest.fixture(scope="module")
def from_hdf_buffer(hdf_file_path, model):
    hdf_buffer = MontecarloRunner.from_hdf(
        hdf_file_path, 'runner', model, Mock_plasma())
    return hdf_buffer


runner_properties = ['nu_bar_estimator',
                     'j_estimator',
                     'last_interaction_in_nu',
                     'last_line_interaction_in_id',
                     'last_line_interaction_out_id',
                     'last_line_interaction_shell_id',
                     'seed', 'inner_boundary_albedo',
                     'enable_reflective_inner_boundary',
                     '_output_nu', '_output_energy']


@pytest.mark.parametrize("attr", runner_properties)
def test_from_hdf_runner(from_hdf_buffer, actual, attr):
    if hasattr(actual, attr):
        npt.assert_almost_equal(getattr(from_hdf_buffer, attr), getattr(
            actual, attr))


runner_quantity_attrs = ['packet_luminosity',
                         'spectrum_frequency', 'sigma_thomson',
                         'montecarlo_virtual_luminosity']


@pytest.mark.parametrize("attr", runner_quantity_attrs)
def test_hdf_runner_quantites(from_hdf_buffer, actual, attr):
    if hasattr(actual, attr):
        assert_quantity_allclose(getattr(from_hdf_buffer, attr), getattr(
            actual, attr))
