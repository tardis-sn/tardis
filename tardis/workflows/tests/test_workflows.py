import pandas as pd
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow
from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow


@pytest.fixture(scope="function")
def v_inner_config(example_configuration_dir):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_config_v_inner.yml"
    )

@pytest.fixture(scope="function")
def run_v_inner(v_inner_config, atomic_data_fname):
    sim = InnerVelocitySolverWorkflow(v_inner_config)
    sim.run()
    return sim

def test_v_inner_solver_workflow(run_v_inner):
    actual = getattr(run_v_inner, 'completed_iterations')
    if hasattr(actual, "value"):
        actual = actual.value
    expected = 10 
    assert actual == expected