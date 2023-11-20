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
def raw_simulation_state(tardis_model_density_config, kurucz_atomic_data):
    return SimulationState.from_config(
        tardis_model_density_config, atom_data=kurucz_atomic_data
    )


@pytest.fixture()
def raw_plasma(
    tardis_model_density_config, raw_simulation_state, kurucz_atomic_data
):
    return assemble_plasma(
        tardis_model_density_config, raw_simulation_state, kurucz_atomic_data
    )


def test_electron_densities(raw_plasma, snapshot_np):
    assert snapshot_np == raw_plasma.electron_densities[8]
    assert snapshot_np == raw_plasma.electron_densities[3]


def test_isotope_number_densities(request, raw_simulation_state, snapshot_np):
    composition = raw_simulation_state.composition
    assert snapshot_np == composition.isotopic_number_density.loc[(28, 56), 0]
    assert snapshot_np == composition.isotopic_number_density.loc[(28, 58), 1]


def test_t_rad(raw_plasma, snapshot_np):
    assert snapshot_np == raw_plasma.t_rad[5]
    assert snapshot_np == raw_plasma.t_rad[3]
