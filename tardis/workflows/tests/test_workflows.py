import pytest
from pathlib import Path
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow
from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow


@pytest.fixture(scope="module")
def config_simulation(example_configuration_dir, atomic_data_fname):
    # this config is intended to match with the one used for simulation_one_loop in the test_simulation.py
    config = Configuration.from_yaml(example_configuration_dir / "tardis_configv1_verysimple.yml")
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)
    return config


@pytest.fixture(scope="module")
def v_inner_config(example_configuration_dir, atomic_data_fname):
    config = Configuration.from_yaml(example_configuration_dir / "tardis_configv1_verysimple.yml")
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    config.model.structure.velocity.start = 5000 * u.km / u.s
    config.model.structure.velocity.num = 50
    config.montecarlo.convergence_strategy["v_inner_boundary"] = {
        "damping_constant": 0.5,
        "threshold": 0.01,
        "type": "damped",
        "store_iteration_properties": True,
    }
    return config


@pytest.fixture(scope="module")
def v_inner_workflow(v_inner_config):
    workflow = InnerVelocitySolverWorkflow(v_inner_config)
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def simple_workflow_one_loop(config_simulation):
    workflow = SimpleTARDISWorkflow(config_simulation)
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def standard_workflow_one_loop(config_simulation):
    workflow = StandardTARDISWorkflow(config_simulation)
    workflow.run()
    return workflow


@pytest.fixture(scope="function")
def simulation_regression_data(regression_data):
    """Fixture to access simulation regression data for test_simulation.py"""

    class SimulationRegressionData:
        def __init__(self, base_regression_data):
            self.base_regression_data = base_regression_data
            self.simulation_regression_dir = (
                base_regression_data.absolute_regression_data_dir.parent
                / "tardis"
                / "simulation"
                / "tests"
                / "test_simulation"
            )

        def get_data(self, attr_root, attr):
            regression_file = self.simulation_regression_dir / f"test_{attr_root}__{attr}__.h5"
            return pd.read_hdf(regression_file)

    return SimulationRegressionData(regression_data)


@pytest.mark.parametrize(
    ["attr_type", "attr"],
    [
        ("plasma_state_iterations", "iterations_w"),
        ("plasma_state_iterations", "iterations_t_rad"),
        ("plasma_state_iterations", "iterations_electron_densities"),
        ("plasma_state_iterations", "iterations_t_inner"),
        # ("plasma_estimates", "nu_bar_estimator"),
        # ("plasma_estimates", "j_estimator"),
        # ("plasma_estimates", "output_nus"),
        # ("plasma_estimates", "output_energies"),
    ],
)
def test_standard_tardis_workflow(
    standard_workflow_one_loop, attr_type, attr, simulation_regression_data
):
    if attr_type == "plasma_estimates":
        if attr in ["nu_bar_estimator", "j_estimator"]:
            actual = getattr(
                standard_workflow_one_loop.transport_state.radfield_mc_estimators,
                attr,
            )
        elif attr in ["output_nus", "output_energies"]:
            actual = getattr(
                standard_workflow_one_loop.transport_state.packet_collection,
                attr,
            )
        else:
            raise ValueError(f"Unknown plasma_estimates attr: {attr}")
        if hasattr(actual, "value"):
            actual = actual.value
        actual = pd.Series(actual)
        expected = simulation_regression_data.get_data(attr_type, attr)
        pd.testing.assert_series_equal(actual, expected, check_exact=False, rtol=1e-3)
    elif attr_type == "plasma_state_iterations":
        attr_data = getattr(standard_workflow_one_loop, attr)
        if hasattr(attr_data, "value"):
            attr_data = attr_data.value
        attr_data = pd.DataFrame(attr_data)
        ref_data = simulation_regression_data.get_data(attr_type, attr)
        pd.testing.assert_frame_equal(attr_data, ref_data, atol=1e-3, rtol=1e-6)
    else:
        raise ValueError(f"Unknown attr_type: {attr_type}")


@pytest.mark.parametrize(
    "attr",
    [
        "t_radiative",
        "dilution_factor",
        "t_inner",
    ],
)
def test_simple_tardis_workflow(simple_workflow_one_loop, standard_workflow_one_loop, attr):
    # this test reply on that the standard workflow pass the test comparing with the regression data
    attr_simple_workflow = getattr(simple_workflow_one_loop.simulation_state, attr)
    attr_standard_workflow = getattr(standard_workflow_one_loop.simulation_state, attr)
    if hasattr(attr_simple_workflow, "value"):
        attr_simple_workflow = attr_simple_workflow.value
    if hasattr(attr_standard_workflow, "value"):
        attr_standard_workflow = attr_standard_workflow.value
    assert np.allclose(attr_simple_workflow, attr_standard_workflow)


@pytest.mark.parametrize(
    "attr",
    [
        "iterations_t_inner",
        "iterations_t_rad",
        "iterations_w",
        "iterations_mean_optical_depth",
        "iterations_v_inner_boundary",
    ],
)
def test_v_inner_solver_workflow(v_inner_workflow, attr, regression_data):
    attr_data = getattr(v_inner_workflow, attr)
    if hasattr(attr_data, "value"):
        attr_data = attr_data.value
    attr_data = pd.DataFrame(attr_data)
    ref_data = regression_data.sync_dataframe(attr_data)
    pd.testing.assert_frame_equal(attr_data, ref_data, atol=1e-3, rtol=1e-6)
