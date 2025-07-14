import pytest
import numpy.testing as npt

from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
from tardis.spectrum.formal_integral.source_function_solver import SourceFunctionSolver

def test_source_function_solver(simulation_verysimple, regression_data, tmp_path):
    """
    Test the output parameters of the SourceFunctionSolver
    """
    sim_state = simulation_verysimple.simulation_state
    plasma = simulation_verysimple.plasma
    transport = simulation_verysimple.transport

    formal_integrator = FormalIntegrator(sim_state, plasma, transport)
    sourceFunction = SourceFunctionSolver(line_interaction_type = formal_integrator.transport.line_interaction_type)
    res = sourceFunction.solve(
        formal_integrator.simulation_state, 
        formal_integrator.opacity_state, 
        formal_integrator.transport, 
        formal_integrator.plasma.atomic_data, 
        formal_integrator.plasma.levels
    )
    att_S_ul, Jred_lu, Jblue_lu, e_dot_u = res[0], res[1], res[2], res[3]


    # check against the regression data
    expected_att_S_ul = regression_data.sync_ndarray(att_S_ul)
    npt.assert_allclose(att_S_ul, expected_att_S_ul)

    expected_Jred_lu = regression_data.sync_ndarray(Jred_lu)
    npt.assert_allclose(Jred_lu, expected_Jred_lu)

    expected_Jblue_lu = regression_data.sync_ndarray(Jblue_lu)
    npt.assert_allclose(Jblue_lu, expected_Jblue_lu)
    
    expected_e_dot_u = regression_data.sync_ndarray(e_dot_u)
    npt.assert_allclose(e_dot_u, expected_e_dot_u)
