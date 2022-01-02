import pytest
import numpy as np
from tardis import constants as c

from copy import deepcopy
import numpy.testing as ntest
from numba import cuda
from numba import njit


import tardis.montecarlo.montecarlo_numba.formal_integral_cuda_test as formal_integral_cuda
import tardis.montecarlo.montecarlo_numba.formal_integral as formal_integral_numba
from tardis.montecarlo.montecarlo_numba.numba_interface import NumbaModel
from tardis.montecarlo.montecarlo_numba.numba_interface import numba_plasma_initialize

from tardis.montecarlo.montecarlo_numba.formal_integral_cuda_functions import cuda_searchsorted_value_right


#All tests need to be changed to call the number version of the function as the expected result

@pytest.mark.parametrize(
    ["nu", "T"],
    [
        (1e14, 1e4),
        (0, 1),
        (1, 1),
    ],
)
def test_intensity_black_body_cuda(nu, T):
    actual = np.zeros(3)
    black_body_caller[1, 3](nu, T, actual)
    
    expected = formal_integral_numba.intensity_black_body(nu, T)
    ntest.assert_almost_equal(actual, expected)
    


@cuda.jit
def black_body_caller(nu, T, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.intensity_black_body_cuda(nu, T)
    


@pytest.mark.parametrize("N", (1e2, 1e3, 1e4, 1e5))
def test_trapezoid_integration_cuda(N):
    actual = np.zeros(4)
    
    h = 1.0
    N = int(N)
    data = np.random.random(N)

    expected = formal_integral_numba.trapezoid_integration(data, h)
    trapezoid_integration_caller[1, 4](data, h, actual)
    
    N_loc = 0
    if N == 1e2:
        N_loc = 0
    elif N == 1e3:
        N_loc = 1
    elif N == 1e4:
        N_loc = 2
    elif N == 1e5:
        N_loc = 3
    
    ntest.assert_almost_equal(actual[N_loc], expected)


@cuda.jit
def trapezoid_integration_caller(data, h, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.trapezoid_integration_cuda(data, h)



@njit(fastmath=True, error_model="numpy", parallel=False)
def calculate_z(r, p, inv_t):
    return np.sqrt(r * r - p * p) * formal_integral_cuda.C_INV * inv_t



TESTDATA_model = [
    {
        "r": np.linspace(1, 2, 3, dtype=np.float64),
    },
    {
        "r": np.linspace(0, 1, 3),
    },
    #{"r": np.linspace(1, 2, 10, dtype=np.float64)},
]
@pytest.fixture(scope="function", params=TESTDATA_model)
def formal_integral_model(request):
    r = request.param["r"]
    model = NumbaModel(
            r[:-1],
            r[1:],
            1/c.c.cgs.value)
    return model



@pytest.mark.parametrize("p", [0.0, 0.5, 1.0])
def test_calculate_z_cuda(formal_integral_model, p):
    actual = np.zeros(3)
    inv_t = 1.0 / formal_integral_model.time_explosion
    size = len(formal_integral_model.r_outer)
    r_outer = formal_integral_model.r_outer 
    for r in r_outer:
        #This is to find which p is being tested, so the assert statement will not crash as
        #actual is an array and needs to be indexed
        p_loc = 0
        if p == 0.0:
            p_loc = 0
        elif p == 0.5:
            p_loc = 1
        elif p == 1.0:
            p_loc = 2
        
        calculate_z_caller[1, 3](r, p, inv_t, actual)
        if p >= r:
            assert actual[p_loc] == 0
        else:
            expected = formal_integral_numba.calculate_z(r, p, inv_t)
            ntest.assert_almost_equal(actual[p_loc], expected)

@cuda.jit
def calculate_z_caller(r, p, inv_t, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.calculate_z_cuda(r, p, inv_t)



@pytest.mark.parametrize("p", [0, 0.5, 1])
def test_populate_z_photosphere(formal_integral_model, p):
    """
    Test the case where p < r[0]
    That means we 'hit' all shells from inside to outside.
    """
    actual = np.zeros(3)
    integrator = formal_integral_cuda.FormalIntegrator(
        formal_integral_model, None, None
    )
    
    size = len(formal_integral_model.r_outer)
    r_inner = formal_integral_model.r_inner
    r_outer = formal_integral_model.r_outer

    p = r_inner[0] * p
    oz = np.zeros_like(r_inner)
    oshell_id = np.zeros_like(oz, dtype=np.int64)
    print(dir(formal_integral_model))
    populate_z_photosphere_caller[1, 3](r_inner, r_outer, formal_integral_model.time_explosion, p, oz, oshell_id, actual)
    #This is to find which p is being tested, so the assert statement will not crash as
    #actual is an array and needs to be indexed
    p_loc = 0
    if p == 0.0:
        p_loc = 0
    elif p == 0.5:
        p_loc = 1
    elif p == 1.0:
        p_loc = 2
    
    expected_oz = np.zeros_like(r_outer)
    expected_oshell_id = np.zeros_like(oshell_id).astype(np.float64)
    formal_integral_numba.populate_z(formal_integral_model, p, expected_oz, expected_oshell_id)
    
    ntest.assert_almost_equal(actual[p_loc], size)

    ntest.assert_allclose(oshell_id, expected_oshell_id)
    
    ntest.assert_allclose(oz, expected_oz, atol=1e-5)

@cuda.jit
def populate_z_photosphere_caller(r_inner, r_outer, time_explosion, p, oz, oshell_id, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.populate_z_cuda(r_inner,
                                                     r_outer,
                                                     time_explosion,
                                                     p, 
                                                     oz, 
                                                     oshell_id)



@pytest.mark.parametrize("p", [1e-5, 1e-3, .1, 0.5, 0.99, 1])
def test_populate_z_shells(formal_integral_model, p):
    """
    Test the case where p > r[0]
    
    oz is redshift
    """
    actual = np.zeros(6)
    
    size = len(formal_integral_model.r_inner)
    r_inner = formal_integral_model.r_inner
    r_outer = formal_integral_model.r_outer
    
    inv_t = 1/formal_integral_model.time_explosion

    p = r_inner[0] + (r_outer[-1] - r_inner[0]) * p
    idx = np.searchsorted(r_outer, p, side="right")

    oz = np.zeros(size * 2)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    offset = size - idx

    expected_N = (offset) * 2
    expected_oz = np.zeros_like(oz)
    expected_oshell_id = np.zeros_like(oshell_id).astype(np.float64)

    # Calculated way to determine which shells get hit
    formal_integral_numba.populate_z(formal_integral_model, p, expected_oz, expected_oshell_id)
    
    populate_z_shells_caller[1,4](r_inner, 
                                  r_outer, 
                                  formal_integral_model.time_explosion, 
                                  p, 
                                  oz, 
                                  oshell_id, 
                                  actual)

    #This is to find which p is being tested, so the assert statement will not crash as
    #actual is an array and needs to be indexed
    
    p_loc = 0
    if p == 1e-5:
        p_loc = 0
    elif p == 0.5:
        p_loc = 1
    elif p == 0.99:
        p_loc = 2
    elif p == 1:
        p_loc = 3

    ntest.assert_almost_equal(actual[p_loc], expected_N)

    ntest.assert_allclose(oshell_id, expected_oshell_id)

    ntest.assert_allclose(oz, expected_oz, atol=1e-5)
    
    ntest.assert_almost_equal(oshell_id, expected_oshell_id)
    
    ntest.assert_almost_equal(oz, expected_oz)

@cuda.jit
def populate_z_shells_caller(r_inner, r_outer, time_explosion, p, oz, oshell_id, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.populate_z_cuda(r_inner,
                                                     r_outer,
                                                     time_explosion,
                                                     p, 
                                                     oz, 
                                                     oshell_id)



@pytest.mark.parametrize(
    "N",
    [
        100,
        1000,
        10000,
    ],
)
def test_calculate_p_values(N):
    #Compare to old calculate p_values
    r = 1.0

    expected = formal_integral_numba.calculate_p_values(r, N)
    
    actual = np.zeros_like(expected, dtype=np.float64)
    actual[::] = formal_integral_cuda.calculate_p_values(r, N)
    
    ntest.assert_allclose(actual, expected)
    
    ntest.assert_almost_equal(actual, expected)



@pytest.mark.parametrize(
    "nu_insert", np.linspace(3e+12, 3e+16, 10)
)
def test_line_search_cuda(nu_insert, verysimple_numba_plasma):
    actual = np.zeros(1)
    expected = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu
   
    expected[0] = formal_integral_numba.line_search(line_list_nu, nu_insert, len(line_list_nu))
    
    line_search_cuda_caller[1, 1](line_list_nu, np.array([nu_insert]), actual)
    
    ntest.assert_almost_equal(actual, expected)

@cuda.jit
def line_search_cuda_caller(line_list_nu, nu_insert, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.line_search_cuda(line_list_nu, nu_insert[x], len(line_list_nu))



@pytest.mark.parametrize(
    "nu_insert", np.linspace(3e+12, 3e+16, 10)
)
def test_reverse_binary_search(nu_insert, verysimple_numba_plasma):
    actual = np.zeros(1)
    expected = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu
    
    if (nu_insert > line_list_nu[0]) or (nu_insert < line_list_nu[len(line_list_nu)-1]):
        raise BoundsError
    
    expected[0] =  len(line_list_nu) - 1 - np.searchsorted(line_list_nu[::-1], nu_insert, side="right")

    cuda_searchsorted_value_right_caller[1, 1](line_list_nu, nu_insert, actual)
    
    ntest.assert_almost_equal(actual, expected)

@cuda.jit
def cuda_searchsorted_value_right_caller(line_list_nu, nu_insert, actual):
    x = cuda.grid(1)
    actual[x] = len(line_list_nu) - 1 - cuda_searchsorted_value_right(line_list_nu[::-1], nu_insert)


    
    
 
from tardis.montecarlo.montecarlo_numba.formal_integral import FormalIntegrator
from tardis.montecarlo.montecarlo_numba.formal_integral_cuda_test import FormalIntegrator as cuda_FormalIntegrator

from tardis.io.atom_data.util import download_atom_data
from tardis import run_tardis
from astropy import units as u
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

from tardis.montecarlo.montecarlo_numba.r_packet import (RPacket)
from tardis.montecarlo.montecarlo_numba.numba_interface import (NumbaModel, 
                                                                NumbaPlasma, 
                                                                numba_plasma_initialize, 
                                                                Estimators)

from basic_types_Copy1 import (get_numba_model,
                         get_r_packet,
                         get_numba_plasma,
                         get_estimator) #This file is in the same test directory for easy access right now. 

from tardis.montecarlo import MontecarloRunner
from tardis.model import Radial1DModel
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.io.util import HDFWriterMixin
from tardis.io.config_reader import ConfigurationError


#This test is a mess right now, I am just trying to see how to start it
@pytest.mark.parametrize(
    ["no_of_packets", "iterations"],
    [
        (40000, 1),
    ],
)
def test_full_formal_integral(no_of_packets, iterations):
    """
    This tests the full formal integral
    """
    download_atom_data('kurucz_cd23_chianti_H_He')
    #This file is a simple file, 1 iteration with only 40000 packets.
    tardis_config = Configuration.from_yaml('tardis/montecarlo/montecarlo_numba/tests/tardis_example_verysimple.yml') 
    #current throws an error here.
    """
    OSError: Atom Data tardis/montecarlo/montecarlo_numba/tests/kurucz_cd23_chianti_H_He.h5 is not found in current path or in TARDIS data repo. tardis/montecarlo/montecarlo_numba/tests/kurucz_cd23_chianti_H_He is also not a standard known TARDIS atom dataset.

tardis/io/atom_data/util.py:48: OSError
    """
    sim1 = Simulation.from_config(tardis_config)
    sim1.run()
    
    basic_estimator1 = get_estimator(sim1)
    basic_estimator2 = get_estimator(sim1)
    
    tau_sobolev_shape1 = basic_estimator1.Edotlu_estimator.shape
    tau_sobolev_shape2 = basic_estimator2.Edotlu_estimator.shape
    
    basic_runner1 = MontecarloRunner.from_config(tardis_config, None, False)
    basic_runner2 = MontecarloRunner.from_config(tardis_config, None, False)

    basic_runner1._initialize_estimator_arrays(sim1.plasma.tau_sobolevs.shape)
    basic_runner2._initialize_estimator_arrays(sim1.plasma.tau_sobolevs.shape)

    basic_runner1._initialize_geometry_arrays(sim1.model)
    basic_runner2._initialize_geometry_arrays(sim1.model)

    basic_runner1._initialize_packets(sim1.model.t_inner.value, 
                                      no_of_packets, 
                                      iterations, 
                                      sim1.model.r_inner[0])
    basic_runner2._initialize_packets(sim1.model.t_inner.value, 
                                      no_of_packets, 
                                      iterations, 
                                      sim1.model.r_inner[0]) 

    basic_runner1.run(sim1.model, sim1.plasma, no_of_packets) #100000 is number of packets
    basic_runner2.run(sim1.model, sim1.plasma, no_of_packets)

    basic_plasma1 = get_numba_plasma(sim1)
    basic_plasma2 = get_numba_plasma(sim1)

    formal_integrator_good = FormalIntegrator(sim1.model, sim1.plasma, basic_runner1)

    formal_integrator_test = cuda_FormalIntegrator(sim1.model, sim1.plasma, basic_runner2)


    res_good = formal_integrator_good.make_source_function()
    att_S_ul_good = res_good[0].flatten(order="F")
    Jred_lu_good = res_good[1].values.flatten(order="F")
    Jblue_lu_good = res_good[2].flatten(order="F")

    res_test = formal_integrator_test.make_source_function()
    att_S_ul_test = res_test[0].flatten(order="F")
    Jred_lu_test = res_test[1].values.flatten(order="F")
    Jblue_lu_test = res_test[2].flatten(order="F")


    formal_integrator_good.generate_numba_objects()

    formal_integrator_test.generate_numba_objects()



    L_test = formal_integrator_test.numba_integrator.formal_integral(
        formal_integrator_test.model.t_inner,
        sim1.runner.spectrum.frequency,
        sim1.runner.spectrum.frequency.shape[0],
        att_S_ul_test,
        Jred_lu_test,
        Jblue_lu_test,
        formal_integrator_test.runner.tau_sobolevs_integ,
        formal_integrator_test.runner.electron_densities_integ,
        formal_integrator_test.points,
    )

    L_good_start = time.monotonic()
    L_good = formal_integrator_good.numba_integrator.formal_integral(
        formal_integrator_good.model.t_inner,
        sim1.runner.spectrum.frequency,
        sim1.runner.spectrum.frequency.shape[0],
        att_S_ul_good,
        Jred_lu_good,
        Jblue_lu_good,
        formal_integrator_good.runner.tau_sobolevs_integ,
        formal_integrator_good.runner.electron_densities_integ,
        formal_integrator_good.points,
    )
    
    ntest.assert_almost_equal(L_test, L_good)













