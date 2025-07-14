import pytest
import numpy.testing as npt

from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver

@pytest.fixture
def source_function_verysimple(simulation_verysimple):
    """
    Generate the source function state from the simulation.
    """
    sim_state = simulation_verysimple.simulation_state
    plasma = simulation_verysimple.plasma
    transport = simulation_verysimple.transport

    formal_integrator = FormalIntegrator(sim_state, plasma, transport)
    sourceFunction = SourceFunctionSolver(line_interaction_type = formal_integrator.transport.line_interaction_type)
    res = sourceFunction.solve(
        formal_integrator.simulation_state, 
        formal_integrator.opacity_state, 
        formal_integrator.transport.transport_state, 
        formal_integrator.plasma.atomic_data, 
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
    expected_e_dot_u = regression_data.sync_ndarray(e_dot_u)
    npt.assert_allclose(e_dot_u, expected_e_dot_u)