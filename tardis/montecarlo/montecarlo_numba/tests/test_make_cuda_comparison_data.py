import pytest
import numpy as np
from tardis import constants as c

from copy import deepcopy
from numba import cuda

import tardis.montecarlo.montecarlo_numba.formal_integral_cuda as formal_integral_cuda
from tardis.montecarlo.montecarlo_numba.numba_interface import NumbaModel

from tardis.montecarlo.montecarlo_numba.formal_integral import FormalIntegrator

# Test cases must also take into account use of a GPU to run. If there is no GPU then the test cases will fail.
GPUs_available = cuda.is_available()

@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(
    ["nu", "T"],
    [
        (1e14, 1e4),
        (0, 1),
        (1, 1),
    ],
)
@pytest.mark.array_compare
def test_intensity_black_body_cuda(nu, T):
    """
    Initializes the test of the cuda version
    against the previous version of itself.
    """
    actual = np.zeros(1)
    black_body_caller[1, 1](nu, T, actual)

    return actual

@cuda.jit
def black_body_caller(nu, T, actual):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.intensity_black_body_cuda(nu, T)



@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(
    "N", (1e2, 1e3, 1e4, 1e5)
)
@pytest.mark.array_compare
def test_trapezoid_integration_cuda(N):
    """
    Initializes the test of the cuda version
    against the previous version of itself.
    """
    actual = np.zeros(1)

    h = 1.0
    N = int(N)
    np.random.seed(12)
    data = np.random.random(N)

    trapezoid_integration_caller[1, 1](data, h, actual)
    return actual

@cuda.jit
def trapezoid_integration_caller(data, h, actual):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.trapezoid_integration_cuda(data, h)


    

TESTDATA_model = [
    {
        "r": np.linspace(1, 2, 3, dtype=np.float64),
    },
    {
        "r": np.linspace(0, 1, 3),
    },
    # {"r": np.linspace(1, 2, 10, dtype=np.float64)},
]


@pytest.fixture(scope="function", params=TESTDATA_model)
def formal_integral_model(request):
    """
    This gets the Numba model to be used in later tests
    """
    r = request.param["r"]
    model = NumbaModel(r[:-1], r[1:], 1 / c.c.cgs.value)
    return model


@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(["p", "p_loc"], [(0.0, 0), (0.5, 1), (1.0, 2)])
@pytest.mark.array_compare
def test_calculate_z_cuda(formal_integral_model, p, p_loc):
    """
    Initializes the test of the cuda version
    against the previous version of itself.
    """
    actual = np.zeros(1)
    inv_t = 1.0 / formal_integral_model.time_explosion
    size = len(formal_integral_model.r_outer)
    r_outer = formal_integral_model.r_outer
    for r in r_outer:
        calculate_z_caller[1, 3](r, p, inv_t, actual)
        return actual

@cuda.jit
def calculate_z_caller(r, p, inv_t, actual):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.calculate_z_cuda(r, p, inv_t)

    

@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(
    ["p", "p_loc"],
    [(1e-5, 0), (1e-3, 1), (0.1, 2), (0.5, 3), (0.99, 4), (1, 5)],
)
@pytest.mark.array_compare
def test_populate_z(formal_integral_model, p, p_loc):
    """
    Initializes the test of the cuda version
    against the previous version of itself
    """
    size = len(formal_integral_model.r_inner)
    oz = np.zeros(size * 2)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    actual = np.zeros(1)
    populate_z_caller[1, 6](
        formal_integral_model.r_inner,
        formal_integral_model.r_outer,
        formal_integral_model.time_explosion,
        p,
        oz,
        oshell_id,
        actual,
    )
    
    big_actual = np.concatenate((actual, oshell_id, oz))
    return big_actual

@cuda.jit
def populate_z_caller(
    r_inner, r_outer, time_explosion, p, oz, oshell_id, actual
):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.populate_z_cuda(
        r_inner, r_outer, time_explosion, p, oz, oshell_id
    )


    
@pytest.mark.parametrize(
    "N",
    [
        100,
        1000,
        10000,
    ],
)
@pytest.mark.array_compare
def test_calculate_p_values(N):
    """
    Initializes the test of the cuda version
    against the previous version of itself.
    """
    r = 1.0
    
    actual = np.zeros(N, dtype=np.float64)
    actual[::] = formal_integral_cuda.calculate_p_values(r, N)

    return actual



@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize("nu_insert", np.linspace(3e12, 3e16, 10))
@pytest.mark.array_compare
def test_line_search_cuda(nu_insert, verysimple_numba_plasma):
    """
    Initializes the test of the cuda version
    against the previous version of itself
    """
    actual = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu

    line_search_cuda_caller[1, 1](line_list_nu, nu_insert, actual)

    return actual

@cuda.jit
def line_search_cuda_caller(line_list_nu, nu_insert, actual):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.line_search_cuda(
        line_list_nu, nu_insert, len(line_list_nu)
    )


    
@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(
    "nu_insert", [*np.linspace(3e12, 3e16, 10), 288786721666522.1]
)
@pytest.mark.array_compare
def test_reverse_binary_search(nu_insert, verysimple_numba_plasma):
    """
    Initializes the test of the cuda version
    against the previous version of itself. 
    The one extra input not included
    in np.linspace a low edge case for testing.
    """
    actual = np.zeros(1)
    #expected = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu

    imin = 0
    imax = len(line_list_nu) - 1

    
    reverse_binary_search_cuda_caller[1, 1](
        line_list_nu, nu_insert, imin, imax, actual
    )

    return actual

@cuda.jit
def reverse_binary_search_cuda_caller(
    line_list_nu, nu_insert, imin, imax, actual
):
    """
    This calls the CUDA function and fills out
    the array
    """
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.reverse_binary_search_cuda(
        line_list_nu, nu_insert, imin, imax
    )


    
# no_of_packets and iterations match what is used by config_verysimple
@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(["no_of_packets", "iterations"], [(200000, 5)])
@pytest.mark.array_compare
def test_full_formal_integral(
    no_of_packets, iterations, config_verysimple, simulation_verysimple
):
    """
    This function initializes the cuda  formal_integrator,
    and the runs them and compares the results to previous 
    versions of itself.
    """
    sim = simulation_verysimple

    formal_integrator_cuda = FormalIntegrator(sim.model, sim.plasma, sim.runner)

    formal_integrator_cuda.interpolate_shells = max(
        2 * formal_integrator_cuda.model.no_of_shells, 80
    )

    res_cuda = formal_integrator_cuda.make_source_function()
    att_S_ul_cuda = res_cuda[0].flatten(order="F")
    Jred_lu_cuda = res_cuda[1].values.flatten(order="F")
    Jblue_lu_cuda = res_cuda[2].flatten(order="F")

    formal_integrator_cuda.generate_numba_objects()

    L_cuda = formal_integrator_cuda.integrator.formal_integral(
        formal_integrator_cuda.model.t_inner,
        sim.runner.spectrum.frequency,
        sim.runner.spectrum.frequency.shape[0],
        att_S_ul_cuda,
        Jred_lu_cuda,
        Jblue_lu_cuda,
        formal_integrator_cuda.runner.tau_sobolevs_integ,
        formal_integrator_cuda.runner.electron_densities_integ,
        formal_integrator_cuda.points,
    )

    return L_cuda
    