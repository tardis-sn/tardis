import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.plasma.standard_plasmas import assemble_plasma


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


def test_electron_densities(raw_plasma, snapshot_np):
    assert snapshot_np == raw_plasma.electron_densities[8]
    assert snapshot_np == raw_plasma.electron_densities[3]


def test_isotope_number_densities(raw_plasma, snapshot_np):
    assert snapshot_np == raw_plasma.isotope_number_density.loc[(28, 56), 0]
    assert snapshot_np == raw_plasma.isotope_number_density.loc[(28, 58), 1]


def test_t_rad(raw_plasma, snapshot_np):
    assert snapshot_np == raw_plasma.t_rad[5]
    assert snapshot_np == raw_plasma.t_rad[3]
