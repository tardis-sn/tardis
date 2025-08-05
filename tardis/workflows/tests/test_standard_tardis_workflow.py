import pandas as pd
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow


@pytest.fixture(scope="module")
def config(example_configuration_dir):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def workflow_one_loop(config, atomic_data_fname):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    workflow = StandardTARDISWorkflow(config)
    workflow.run()
    return workflow




@pytest.mark.parametrize(
    "attr",
    [
        "iterations_w",
        "iterations_t_rad",
        "iterations_electron_densities",
        "iterations_t_inner",
    ],
)
def test_plasma_state_iterations(workflow_one_loop, attr, regression_data):
    actual = getattr(workflow_one_loop, attr)
    if hasattr(actual, "value"):
        actual = actual.value
    actual = pd.DataFrame(actual)
    expected = regression_data.sync_dataframe(actual)
    pd.testing.assert_frame_equal(actual, expected, rtol=1e-5, atol=1e-8)


@pytest.mark.parametrize(
    "attr",
    [
        "nu_bar_estimator",
        "j_estimator",
        "t_radiative",
        "dilution_factor",
        "output_nus",
        "output_energies",
    ],
)
def test_plasma_estimates(workflow_one_loop, attr, regression_data):
    if attr in ["nu_bar_estimator", "j_estimator"]:
        actual = getattr(
            workflow_one_loop.transport_state.radfield_mc_estimators,
            attr,
        )
    elif attr in ["t_radiative", "dilution_factor"]:
        actual = getattr(workflow_one_loop.simulation_state, attr)
    elif attr in ["output_nus", "output_energies"]:
        actual = getattr(
            workflow_one_loop.transport_state.packet_collection,
            attr,
        )
    else:
        actual = getattr(workflow_one_loop.transport, attr)

    if hasattr(actual, "value"):
        actual = actual.value
    actual = pd.Series(actual)
    expected = regression_data.sync_dataframe(actual)
    pd.testing.assert_series_equal(actual, expected, rtol=1e-5, atol=1e-8)