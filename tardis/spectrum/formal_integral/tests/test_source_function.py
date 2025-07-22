import pytest
import numpy.testing as npt
from copy import deepcopy

from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.simulation import Simulation
from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver

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

    formal_integrator = FormalIntegrator(sim_state, plasma, transport)
    sourceFunction = SourceFunctionSolver(formal_integrator.transport.line_interaction_type, 
                                          formal_integrator.atomic_data)
    res = sourceFunction.solve(
        formal_integrator.simulation_state, 
        formal_integrator.opacity_state, 
        formal_integrator.transport.transport_state,
        formal_integrator.plasma.levels
    )
    return res


def test_att_S_ul(source_function_verysimple, request):
    """
    Test the attenuated source function
    """
    regression_data = RegressionData(request)

    att_S_ul = source_function_verysimple.att_S_ul
    expected_att_S_ul = regression_data.sync_ndarray(att_S_ul)
    npt.assert_allclose(att_S_ul, expected_att_S_ul)

def test_Jred_lu(source_function_verysimple, request):
    """
    Test the red end line estimator
    """
    regression_data = RegressionData(request)

    Jred_lu = source_function_verysimple.Jred_lu
    expected_Jred_lu = regression_data.sync_ndarray(Jred_lu)
    npt.assert_allclose(Jred_lu, expected_Jred_lu)

def test_Jblue_lu(source_function_verysimple, request):
    """
    Test the blue end line estimator
    """
    regression_data = RegressionData(request)
    
    Jblue_lu = source_function_verysimple.Jblue_lu
    expected_Jblue_lu = regression_data.sync_ndarray(Jblue_lu)
    npt.assert_allclose(Jblue_lu, expected_Jblue_lu)

def test_e_dot_u(source_function_verysimple, request):
    """
    Test the energy density estimator
    """
    regression_data = RegressionData(request)
    
    e_dot_u = source_function_verysimple.e_dot_u
    expected_e_dot_u = regression_data.sync_dataframe(e_dot_u)
    npt.assert_allclose(e_dot_u, expected_e_dot_u)