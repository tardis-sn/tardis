import numpy as np
import pandas as pd
import pytest
from astropy import units as u

from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow
from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow


@pytest.fixture(scope="module")
def v_inner_config(config_verysimple_for_simulation_one_loop):
    config_verysimple_for_simulation_one_loop.model.structure.velocity.start = (
        5000 * u.km / u.s
    )
    config_verysimple_for_simulation_one_loop.model.structure.velocity.num = 50
    config_verysimple_for_simulation_one_loop.montecarlo.convergence_strategy[
        "v_inner_boundary"
    ] = {
        "damping_constant": 0.5,
        "threshold": 0.01,
        "type": "damped",
        "store_iteration_properties": True,
    }
    return config_verysimple_for_simulation_one_loop


@pytest.fixture(scope="module")
def v_inner_workflow(v_inner_config):
    workflow = InnerVelocitySolverWorkflow(v_inner_config)
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def simple_workflow_one_loop(config_verysimple_for_simulation_one_loop):
    workflow = SimpleTARDISWorkflow(config_verysimple_for_simulation_one_loop)
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def standard_workflow_one_loop(config_verysimple_for_simulation_one_loop):
    workflow = StandardTARDISWorkflow(config_verysimple_for_simulation_one_loop)
    workflow.run()
    return workflow


@pytest.mark.parametrize(
    ["attr_type", "attr"],
    [
        ("plasma_state_iterations", "iterations_w"),
        ("plasma_state_iterations", "iterations_t_rad"),
        ("plasma_state_iterations", "iterations_electron_densities"),
        ("plasma_state_iterations", "iterations_t_inner"),
        ("plasma_estimates", "nu_bar_estimator"),
        ("plasma_estimates", "j_estimator"),
        ("plasma_estimates", "output_nus"),
        ("plasma_estimates", "output_energies"),
    ],
)
def test_standard_tardis_workflow_against_run_tardis(
    standard_workflow_one_loop, attr_type, attr, regression_data
):
    ref_file = (
        regression_data.regression_data_path
        / "tardis"
        / "simulation"
        / "tests"
        / "test_simulation"
        / f"test_{attr_type}__{attr}__.h5"
    )
    ref_data = pd.read_hdf(ref_file)
    if attr_type == "plasma_estimates":
        if attr in ["nu_bar_estimator", "j_estimator"]:
            # Map old attribute names to new ones
            attr_map = {
                "nu_bar_estimator": "mean_frequency",
                "j_estimator": "mean_intensity_total",
            }
            attr_data = getattr(
                standard_workflow_one_loop.transport_state.estimators_bulk,
                attr_map[attr],
            )
        elif attr in ["output_nus", "output_energies"]:
            attr_data = getattr(
                standard_workflow_one_loop.transport_state.packet_collection,
                attr,
            )
        else:
            raise ValueError(f"Unknown plasma_estimates attr: {attr}")
        if hasattr(attr_data, "value"):
            attr_data = attr_data.value
        attr_data = pd.Series(attr_data)
        pd.testing.assert_series_equal(
            attr_data, ref_data, check_exact=False, rtol=1e-6
        )
    elif attr_type == "plasma_state_iterations":
        attr_data = getattr(standard_workflow_one_loop, attr)
        if hasattr(attr_data, "value"):
            attr_data = attr_data.value
        attr_data = pd.DataFrame(attr_data)
        pd.testing.assert_frame_equal(attr_data, ref_data, atol=0, rtol=1e-14)
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
def test_simple_tardis_workflow_against_standard_workflow(
    simple_workflow_one_loop, standard_workflow_one_loop, attr
):
    # this test cross compare the simple workflow results against the standard workflow results for now,
    # since by default the simple workflow does not save the iteration properties
    attr_simple_workflow = getattr(
        simple_workflow_one_loop.simulation_state, attr
    )
    attr_standard_workflow = getattr(
        standard_workflow_one_loop.simulation_state, attr
    )
    if hasattr(attr_simple_workflow, "value"):
        attr_simple_workflow = attr_simple_workflow.value
    if hasattr(attr_standard_workflow, "value"):
        attr_standard_workflow = attr_standard_workflow.value
    assert np.allclose(
        attr_simple_workflow, attr_standard_workflow, atol=0, rtol=1e-14
    )


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
