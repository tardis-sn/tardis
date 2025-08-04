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


@pytest.fixture(scope="module")
def simulation_regression_data(tardis_regression_path):
    regression_data_path = tardis_regression_path / "tardis" / "simulation" / "tests" / "test_simulation"
    h5_file_path = regression_data_path / "test_plasma_state_iterations__iterations_w__.h5"
    
    store = pd.HDFStore(h5_file_path, mode="r")
    yield store
    store.close()


def test_plasma_state_iterations_w(workflow_one_loop, simulation_regression_data):
    actual = workflow_one_loop.iterations_w
    if hasattr(actual, "value"):
        actual = actual.value
    actual = pd.DataFrame(actual)
    
    expected = pd.read_hdf(simulation_regression_data, key="data")
    
    pd.testing.assert_frame_equal(actual, expected, rtol=1e-5, atol=1e-8)