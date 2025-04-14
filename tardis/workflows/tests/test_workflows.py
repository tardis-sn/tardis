import pytest
import numpy as np
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow
from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow

@pytest.fixture(scope="module")
def simulation_simple_workflow(config_verysimple, atomic_data_fname):
    sim = SimpleTARDISWorkflow(config_verysimple, csvy=True)
    sim.run()
    return sim

@pytest.fixture(scope="module")
def simulation_standard_workflow(config_verysimple, atomic_data_fname):
    sim = StandardTARDISWorkflow(config_verysimple, csvy=True)
    sim.run()
    return sim

@pytest.fixture(scope="module")
def v_inner_config(example_configuration_dir):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_config_v_inner.yml"
    )

@pytest.fixture(scope="module")
def simulation_v_inner_workflow(v_inner_config, atomic_data_fname):
    sim = InnerVelocitySolverWorkflow(v_inner_config)
    sim.run()
    return sim

def test_v_inner_solver_workflow_iterations(simulation_v_inner_workflow):
    assert hasattr(simulation_v_inner_workflow, "completed_iterations")
    assert simulation_v_inner_workflow.completed_iterations == 10

@pytest.mark.parametrize(["key","expected"],[("t_inner",10684.13082565328 * u.K),("v_inner",1970000000.0 * u.cm / u.s)])
def test_v_inner_solver_workflow_convergence(key, expected, simulation_v_inner_workflow):
    actual = getattr(simulation_v_inner_workflow.simulation_state, key)
    if not actual.isscalar:
        actual = actual[-1]    
    assert np.allclose(actual.to_value(), expected.to_value())

@pytest.mark.parametrize("key",["t_radiative","dilution_factor","t_inner"])
def test_simple_tardis_workflow(key, simulation_simple_workflow, simulation_verysimple):
    actual = getattr(simulation_simple_workflow.simulation_state, key)
    actual = actual.to_value() if hasattr(actual, "to_value") else actual
    expected = getattr(simulation_verysimple.simulation_state, key)
    expected = expected.to_value() if hasattr(expected, "to_value") else expected
    assert np.allclose(actual, expected)

@pytest.mark.parametrize("key",["t_radiative","dilution_factor","t_inner"])
def test_standard_tardis_workflow(key, simulation_standard_workflow, simulation_verysimple):
    actual = getattr(simulation_standard_workflow.simulation_state, key)
    actual = actual.to_value() if hasattr(actual, "to_value") else actual
    expected = getattr(simulation_verysimple.simulation_state, key)
    expected = expected.to_value() if hasattr(expected, "to_value") else expected
    assert np.allclose(actual, expected)