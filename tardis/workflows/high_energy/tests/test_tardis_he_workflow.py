"""
Test suite for TARDIS High Energy (HE) Workflow

This module contains comprehensive tests for the TARDISHEWorkflow class,
following TARDIS testing conventions and using regression data.
"""

import pytest
from numpy.testing import assert_allclose

from tardis.workflows.high_energy.tardis_he_workflow import TARDISHEWorkflow


@pytest.fixture(scope="session")
def he_workflow_instance_session(atomic_dataset, he_test_config):
    """Create a TARDISHEWorkflow instance for testing (session-scoped)."""
    return TARDISHEWorkflow(atomic_dataset, he_test_config, config_type="yaml")


@pytest.fixture(scope="session")
def he_workflow_minimal_run_params():
    """Minimal parameters for HE workflow run."""
    return {
        "time_start": 0.1,
        "time_end": 10.0,
        "number_of_packets": 1000,
        "time_steps": 5,
        "time_space": "log",
        "seed": 23,
        "fp": 0.5,
        "spectrum_bins": 100,
        "grey_opacity": -1,
        "legacy": False,
        "legacy_atom_data": None,
    }


@pytest.fixture(scope="session")
def he_workflow_result(he_workflow_instance_session, he_workflow_minimal_run_params):
    """Run the HE workflow once and cache the result for multiple tests."""
    return he_workflow_instance_session.run(**he_workflow_minimal_run_params)


@pytest.fixture(scope="session")
def he_workflow_time_spacing_params():
    """Parameter sets for time spacing tests."""
    return {
        "linear": {
            "time_start": 1.0,
            "time_end": 20.0,
            "number_of_packets": 500,
            "time_steps": 4,
            "time_space": "linear",
            "seed": 23,
            "fp": 0.5,
            "spectrum_bins": 50,
            "grey_opacity": -1,
            "legacy": False,
            "legacy_atom_data": None,
        },
        "log": {
            "time_start": 1.0,
            "time_end": 20.0,
            "number_of_packets": 500,
            "time_steps": 4,
            "time_space": "log",
            "seed": 23,
            "fp": 0.5,
            "spectrum_bins": 50,
            "grey_opacity": -1,
            "legacy": False,
            "legacy_atom_data": None,
        },
    }


@pytest.fixture(scope="session")
def he_workflow_time_spacing_results(
    he_workflow_instance_session, he_workflow_time_spacing_params
):
    """Run the HE workflow with different time spacing parameters and cache the results."""
    results = {}
    for time_space, params in he_workflow_time_spacing_params.items():
        results[time_space] = he_workflow_instance_session.run(**params)
    return results


@pytest.mark.parametrize(
    "output_he",
    [
        "escape_energy",
        "escape_energy_cosi",
        "packets_escaped",
        "gamma_ray_deposited_energy",
        "total_deposited_energy",
        "positron_energy",
    ],
)
def test_he_workflow_all_outputs_regression(
    he_workflow_result, regression_data, output_he
):
    """
    Test each main HE workflow output (DataFrame) against regression data.
    """
    output = getattr(he_workflow_result, output_he)
    expected = regression_data.sync_dataframe(output, key=output_he)
    assert_allclose(output.values, expected.values, rtol=1e-5)
