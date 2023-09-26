import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.plasma.standard_plasmas import assemble_plasma
from numpy.testing import assert_almost_equal


@pytest.fixture
def tardis_model_density_config(example_model_file_dir):
    return Configuration.from_yaml(
        example_model_file_dir / "tardis_configv1_tardis_model_format.yml"
    )


@pytest.fixture()
def raw_model(tardis_model_density_config):
    return SimulationState.from_config(tardis_model_density_config)


@pytest.fixture()
def raw_plasma(tardis_model_density_config, raw_model, kurucz_atomic_data):
    return assemble_plasma(
        tardis_model_density_config, raw_model, kurucz_atomic_data
    )


def test_electron_densities(raw_plasma):
    assert_almost_equal(raw_plasma.electron_densities[8], 2.72e14)
    assert_almost_equal(raw_plasma.electron_densities[3], 2.6e14)


def test_isotope_number_densities(raw_plasma):
    assert_almost_equal(
        raw_plasma.isotope_number_density.loc[(28, 56), 0], 9688803936.317898
    )
    assert_almost_equal(
        raw_plasma.isotope_number_density.loc[(28, 58), 1], 13097656958.746628
    )


def test_t_rad(raw_plasma):
    assert_almost_equal(raw_plasma.t_rad[5], 76459.592)
    assert_almost_equal(raw_plasma.t_rad[3], 76399.042)
