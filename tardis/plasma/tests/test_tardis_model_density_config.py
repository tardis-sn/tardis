import pytest
import os
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma
from numpy.testing import assert_almost_equal

data_path = os.path.join("tardis", "io", "tests", "data")


@pytest.fixture
def tardis_model_density_config():
    filename = "tardis_configv1_tardis_model_format.yml"
    return Configuration.from_yaml(os.path.join(data_path, filename))


@pytest.fixture()
def raw_model(tardis_model_density_config):
    return Radial1DModel.from_config(tardis_model_density_config)


@pytest.fixture()
def raw_plasma(tardis_model_density_config, raw_model, kurucz_atomic_data):
    return assemble_plasma(
        tardis_model_density_config, raw_model, kurucz_atomic_data
    )


def test_electron_densities(raw_plasma):
    assert_almost_equal(raw_plasma.electron_densities[8], 2.72e14)
    assert_almost_equal(raw_plasma.electron_densities[3], 2.6e14)


def test_t_rad(raw_plasma):
    assert_almost_equal(raw_plasma.t_rad[5], 76459.592)
    assert_almost_equal(raw_plasma.t_rad[3], 76399.042)
