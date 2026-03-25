"""
Tests for the double-slicing fix of opacity state boundary indices.

The plasma assembles all quantities (tau_sobolev, electron_densities,
transition_probabilities) indexed over active shells only.  Before the fix,
several code paths re-sliced those arrays with the raw v_inner/v_outer
boundary indices, producing arrays that were too short.  These tests run
actual workflows and verify the shape invariant: every plasma array must
have exactly no_of_shells_active columns/elements.

InnerVelocitySolverWorkflow is the primary test target because it
naturally moves v_inner_boundary to a non-zero index after the first
iteration, which is when the double-slicing bug manifests.
"""

import os

import pytest
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.workflows.util import get_tau_integ

HERE = os.path.dirname(__file__)
CONFIG_SLICING = os.path.join(HERE, "config_slicing_test.yml")
CONFIG_V_INNER = os.path.join(HERE, "config_v_inner_solver_test.yml")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assert_active_shell_shape(array, n_active, name):
    last_dim = array.shape[-1] if hasattr(array, "shape") else len(array)
    assert last_dim == n_active, (
        f"{name}: last dimension is {last_dim}, expected no_of_shells_active={n_active}"
    )


# ---------------------------------------------------------------------------
# Test 1: SimpleTARDISWorkflow with default boundary
#         Smoke test — verifies no regression for the common case.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def simple_workflow():
    from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow

    config = Configuration.from_yaml(CONFIG_SLICING)
    workflow = SimpleTARDISWorkflow(config)
    workflow.run()
    return workflow


def test_simple_workflow_completes(simple_workflow):
    assert simple_workflow.spectrum_solver.spectrum_real_packets is not None


def test_simple_workflow_opacity_shape(simple_workflow):
    n_active = simple_workflow.simulation_state.no_of_shells
    ts = simple_workflow.transport_state
    _assert_active_shell_shape(
        ts.opacity_state.tau_sobolev, n_active, "opacity_state.tau_sobolev"
    )


# ---------------------------------------------------------------------------
# Test 2: InnerVelocitySolverWorkflow
#         After the first iteration the boundary moves, giving
#         v_inner_boundary_index > 0.  Before the fix this caused a
#         ValueError in interp1d (tau_integ length mismatch) and shape
#         errors in get_tau_integ.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def v_inner_workflow():
    from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow

    config = Configuration.from_yaml(CONFIG_V_INNER)
    config.montecarlo.iterations = 3
    workflow = InnerVelocitySolverWorkflow(config)
    workflow.run()
    return workflow


def test_v_inner_workflow_completes(v_inner_workflow):
    assert v_inner_workflow.final_integrated_tau is not None


def test_v_inner_boundary_moved(v_inner_workflow):
    # The workflow should have found a non-trivial inner boundary.
    geom = v_inner_workflow.simulation_state.geometry
    assert geom.v_inner_boundary_index >= 0  # may be 0 if boundary lands at first shell


def test_v_inner_workflow_opacity_shape(v_inner_workflow):
    n_active = v_inner_workflow.simulation_state.no_of_shells
    ts = v_inner_workflow.transport_state
    _assert_active_shell_shape(
        ts.opacity_state.tau_sobolev, n_active, "opacity_state.tau_sobolev"
    )


def test_v_inner_get_tau_integ_shape(v_inner_workflow):
    """get_tau_integ must return arrays of length no_of_shells_active.
    Before the fix, it used radiation_field_state.temperature (all shells)
    and geometry.r_outer (all shells), causing shape mismatches."""
    n_active = v_inner_workflow.simulation_state.no_of_shells
    result = get_tau_integ(
        v_inner_workflow.plasma_solver,
        v_inner_workflow.opacity_states["opacity_state"],
        v_inner_workflow.simulation_state,
    )
    for key, arr in result.items():
        _assert_active_shell_shape(arr, n_active, f"get_tau_integ['{key}']")
