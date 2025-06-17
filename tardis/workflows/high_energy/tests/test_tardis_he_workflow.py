import pytest
from numpy.testing import assert_allclose

from tardis.workflows.high_energy.tardis_he_workflow import TARDISHEWorkflow


@pytest.fixture(scope="session")
def he_workflow_instance_session(atomic_dataset, he_test_config):
    """Create a TARDISHEWorkflow instance for testing (session-scoped)."""
    return TARDISHEWorkflow(atomic_dataset, he_test_config, config_type="csvy")


@pytest.fixture(scope="session")
def he_workflow_minimal_run_params(atomic_dataset):
    """Minimal parameters for HE workflow run."""
    return {
        "time_start": 0.1,
        "time_end": 100.0,
        "number_of_packets": int(3e5),
        "time_steps": 20,
        "time_space": "log",
        "seed": 1993,
        "fp": 1.9,
        "spectrum_bins": 1000,
        "grey_opacity": -1,
        "legacy": True,
        "legacy_atom_data": atomic_dataset,
    }


@pytest.fixture(scope="session")
def he_workflow_result(he_workflow_instance_session, he_workflow_minimal_run_params):
    """Run the HE workflow once and cache the result for multiple tests."""
    return he_workflow_instance_session.run(**he_workflow_minimal_run_params)


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
