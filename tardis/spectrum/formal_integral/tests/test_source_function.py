from copy import deepcopy

import numpy.testing as npt
import pytest

from tardis.simulation import Simulation
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver

SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL = 1e-13

config_line_modes = ["downbranch", "macroatom"]


@pytest.fixture(scope="module", params=config_line_modes)
def source_function_verysimple(request, config_verysimple, atomic_dataset):
    """
    Generate the source function state from the simulation.
    """

    config_verysimple["plasma"]["line_interaction_type"] = request.param
    print(f"Using line interaction type: {request.param}")
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.last_no_of_packets = 4000
    sim.run_final()

    sim_state = sim.simulation_state
    plasma = sim.plasma
    transport = sim.transport

    integrator_settings = sim.spectrum_solver.integrator_settings
    formal_integrator = FormalIntegralSolver(
        integrator_settings.points,
        integrator_settings.interpolate_shells,
        getattr(integrator_settings, "method", None),
    )
    opacity_state = formal_integrator.setup(
        transport, sim.opacity_state, sim.macro_atom_state
    )
    atomic_data = plasma.atomic_data
    source_function_solver = SourceFunctionSolver(
        transport.line_interaction_type
    )
    source_function_state = source_function_solver.solve(
        sim_state,
        opacity_state,
        transport.transport_state,
        atomic_data,
        sim.macro_atom_state,
    )
    return source_function_state


def test_att_S_ul(source_function_verysimple, regression_data):
    """
    Test att_S_ul
    """
    att_S_ul = source_function_verysimple.att_S_ul
    expected_att_S_ul = regression_data.sync_ndarray(att_S_ul)
    npt.assert_allclose(
        att_S_ul.mean(axis=0),
        expected_att_S_ul.mean(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )
    npt.assert_allclose(
        att_S_ul.std(axis=0),
        expected_att_S_ul.std(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )


def test_Jred_lu(source_function_verysimple, regression_data):
    """
    Test Jred_lu
    """
    Jred_lu = source_function_verysimple.Jred_lu
    expected_Jred_lu = regression_data.sync_ndarray(Jred_lu)
    npt.assert_allclose(
        Jred_lu.mean(axis=0),
        expected_Jred_lu.mean(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )
    npt.assert_allclose(
        Jred_lu.std(axis=0),
        expected_Jred_lu.std(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )


def test_Jblue_lu(source_function_verysimple, regression_data):
    """
    Test Jblue_lu
    """
    Jblue_lu = source_function_verysimple.Jblue_lu
    expected_Jblue_lu = regression_data.sync_ndarray(Jblue_lu)
    npt.assert_allclose(
        Jblue_lu.mean(axis=0),
        expected_Jblue_lu.mean(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )
    npt.assert_allclose(
        Jblue_lu.std(axis=0),
        expected_Jblue_lu.std(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )


def test_e_dot_u(source_function_verysimple, regression_data):
    """
    Test e_dot_u
    """
    e_dot_u = source_function_verysimple.e_dot_u
    expected_e_dot_u = regression_data.sync_dataframe(e_dot_u)
    npt.assert_allclose(
        e_dot_u.mean(axis=0),
        expected_e_dot_u.mean(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )
    npt.assert_allclose(
        e_dot_u.std(axis=0),
        expected_e_dot_u.std(axis=0),
        rtol=SOURCE_FUNCTION_FORMAL_INTEGRAL_RTOL,
        atol=0,
    )
