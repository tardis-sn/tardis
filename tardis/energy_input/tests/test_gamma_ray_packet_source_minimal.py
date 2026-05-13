"""Test the gamma ray packet source using actual workflow configurations."""

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

import tardis
from tardis.energy_input.gamma_ray_channel import (
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_decay_df,
    create_isotope_dicts,
    time_evolve_cumulative_decay,
)
from tardis.transport.montecarlo.packet_source.high_energy import GammaRayPacketSource
from tardis.energy_input.main_gamma_ray_loop import get_effective_time_array
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState


@pytest.fixture(scope="module")
def simulation_state(atomic_dataset):
    """Create a simulation state from the HE workflow config."""

    # Get tardis root directory
    tardis_root = Path(tardis.__path__[0])
    config_path = (
        tardis_root
        / "workflows"
        / "high_energy"
        / "tests"
        / "data"
        / "tardis_he_test_config.yml"
    )

    return SimulationState.from_csvy(Configuration.from_yaml(str(config_path)))


@pytest.fixture(scope="module")
def total_decays_df(simulation_state):
    """
    Create total decays dataframe using the same method as the workflow.
    """
    isotopic_mass_fraction = simulation_state.composition.nuclide_mass_fraction
    cell_masses = (simulation_state.volume * simulation_state.density).to(u.g)

    isotope_dicts = create_isotope_dicts(isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dicts)

    total_decays = calculate_total_decays(inventories_dict, 86400.0 * u.s)

    return total_decays


@pytest.fixture(scope="module")
def isotope_decay_df(total_decays_df, atomic_dataset):
    """
    Create isotope decay dataframe from total decays and gamma ray lines.
    """
    gamma_ray_lines = atomic_dataset.decay_radiation_data
    isotope_decay_df = create_isotope_decay_df(total_decays_df, gamma_ray_lines)

    return isotope_decay_df


@pytest.fixture(scope="module")
def effective_time_array():
    """Create effective time array for testing."""
    times, effective_times = get_effective_time_array(0.1, 1.0, "linear", 2)
    return times, effective_times


@pytest.fixture(scope="module")
def cumulative_decays_df(simulation_state, atomic_dataset, effective_time_array):
    """
    Create cumulative decays dataframe that evolves over time,
    similar to workflow's time_evolve_cumulative_decay_expanded.
    """
    isotopic_mass_fraction = simulation_state.composition.nuclide_mass_fraction
    cell_masses = (simulation_state.volume * simulation_state.density).to(u.g)
    gamma_ray_lines = atomic_dataset.decay_radiation_data

    times, _ = effective_time_array

    cumulative_decays_df = time_evolve_cumulative_decay(
        isotopic_mass_fraction, cell_masses, gamma_ray_lines, times
    )

    return cumulative_decays_df


@pytest.fixture(scope="module")
def packet_source_params(
    simulation_state, isotope_decay_df, effective_time_array, cumulative_decays_df
):
    """Calculate parameters for packet source from simulation state."""
    inner_velocities = simulation_state.v_inner.to("cm/s").value
    outer_velocities = simulation_state.v_outer.to("cm/s").value

    times, effective_times = effective_time_array

    return {
        "cumulative_decays_df": cumulative_decays_df,
        "isotope_decay_df": isotope_decay_df,
        "positronium_fraction": 0.0,
        "inner_velocities": inner_velocities,
        "outer_velocities": outer_velocities,
        "times": times,
        "effective_times": effective_times,
        "base_seed": 1963,
    }


@pytest.fixture(scope="module")
def gamma_ray_source(packet_source_params):
    """Create a new GammaRayPacketSource using workflow data."""
    return GammaRayPacketSource(**packet_source_params)


@pytest.fixture(scope="module")
def legacy_energy_per_packet(gamma_ray_source):
    """Calculate legacy energy per packet from isotope_decay_df."""
    gamma_df = gamma_ray_source.isotope_decay_df[
        gamma_ray_source.isotope_decay_df["radiation"] == "g"
    ]
    total_energy_gamma = gamma_df["decay_energy_erg"].sum()
    n_packets = 300000  # Number of packets to create

    if total_energy_gamma == 0:
        pytest.skip(
            "No gamma ray decays in the time interval, skipping legacy energy calculation."
        )

    return total_energy_gamma / n_packets


@pytest.fixture(scope="module")
def created_packets_data_legacy(
    gamma_ray_source, cumulative_decays_df, legacy_energy_per_packet
):
    """Creates gamma ray packets using legacy energy calculation from isotope_decay_df."""
    n_packets = 300000  # Number of packets to create
    packets = gamma_ray_source.create_packets(
        cumulative_decays_df,
        number_of_packets=n_packets,
        legacy_energy_per_packet=legacy_energy_per_packet,
    )
    return packets, n_packets


@pytest.fixture(scope="module")
def created_packets_data_cumulative(gamma_ray_source, cumulative_decays_df):
    """Creates gamma ray packets using energy calculation from cumulative_decays_df."""
    # Check if there's any gamma ray energy to create packets from
    gamma_df = gamma_ray_source.cumulative_decays_df[
        gamma_ray_source.cumulative_decays_df["radiation"] == "g"
    ]
    total_energy_gamma = gamma_df["decay_energy_erg"].sum()

    if total_energy_gamma == 0:
        pytest.skip(
            "No gamma ray decays in the time interval, skipping packet creation dependent tests."
        )

    n_packets = 300000  # Number of packets to create
    packets = gamma_ray_source.create_packets(
        cumulative_decays_df, number_of_packets=n_packets
    )
    return packets, n_packets


# Attributes to test for GXPacketCollection, with their expected type format
# (name, type_string where '[:]' indicates an array)
GX_COLLECTION_ATTR_TEST = [
    "location",
    "direction",
    "energy_rf",
    "energy_cmf",
    "nu_rf",
    "nu_cmf",
    "status",
    "shell",
    "time_start",
    "time_index",
]


@pytest.mark.parametrize("gx_collection_attr", GX_COLLECTION_ATTR_TEST)
def test_gamma_ray_packet_properties_legacy_energy(
    gx_collection_attr, created_packets_data_legacy, regression_data
):
    """Test properties of the GXPacketCollection using legacy energy calculation from isotope_decay_df."""
    packets, _n_packets = created_packets_data_legacy

    actual_value = getattr(packets, gx_collection_attr)

    # The distinction based on type_str for array vs scalar is removed as per instruction.
    # We now pass actual_value directly to sync_ndarray.
    value_to_sync = actual_value

    # The key argument is removed as per instruction.
    expected_synced_value = regression_data.sync_ndarray(value_to_sync)

    assert_allclose(value_to_sync, expected_synced_value, rtol=1e-7)


@pytest.mark.parametrize("gx_collection_attr", GX_COLLECTION_ATTR_TEST)
def test_gamma_ray_packet_properties_cumulative_energy(
    gx_collection_attr, created_packets_data_cumulative, regression_data
):
    """Test properties of the GXPacketCollection using energy calculation from cumulative_decays_df."""
    packets, _n_packets = created_packets_data_cumulative

    actual_value = getattr(packets, gx_collection_attr)

    # The distinction based on type_str for array vs scalar is removed as per instruction.
    # We now pass actual_value directly to sync_ndarray.
    value_to_sync = actual_value

    # The key argument is removed as per instruction.
    expected_synced_value = regression_data.sync_ndarray(value_to_sync)

    assert_allclose(value_to_sync, expected_synced_value, rtol=1e-7)
