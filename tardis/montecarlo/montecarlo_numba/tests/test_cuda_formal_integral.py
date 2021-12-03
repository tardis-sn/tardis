import pytest
import numpy as np
from tardis import constants as c

from copy import deepcopy
import numpy.testing as ntest
from numba import cuda

from tardis.util.base import intensity_black_body
import tardis.montecarlo.montecarlo_numba.formal_integral_cuda_test as formal_integral_cuda
from tardis.montecarlo.montecarlo_numba.numba_interface import NumbaModel

from tardis.montecarlo.montecarlo_numba.formal_integral_cuda_functions import cuda_searchsorted_value_right

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
    
    expected = intensity_black_body(nu, T)
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

    #actual = func(data, h)
    expected = np.trapz(data)
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
    

def calculate_z(r, p, inv_t):
    return np.sqrt(r * r - p * p) * formal_integral_cuda.C_INV * inv_t

TESTDATA = [
    {
        "r": np.linspace(1, 2, 3, dtype=np.float64),
    },
    {
        "r": np.linspace(0, 1, 3),
    },
    #{"r": np.linspace(1, 2, 10, dtype=np.float64)},
]


@pytest.fixture(scope="function", params=TESTDATA)
def formal_integral_model(request):
    r = request.param["r"]
    model = NumbaModel(
            r[:-1],
            r[1:],
            1/c.c.cgs.value)
    return model


@pytest.mark.parametrize(
        'p', [0.0, 0.5, 1.0]
)
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
            desired = calculate_z(r, p, inv_t)
            ntest.assert_almost_equal(actual[p_loc], desired)

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
    #func = formal_integral_cuda.populate_z_cuda
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
    
    #N = func(formal_integral_model, p, oz, oshell_id)
    ntest.assert_almost_equal(actual[p_loc], size)

    ntest.assert_allclose(oshell_id, np.arange(0, size, 1))

    ntest.assert_allclose(oz, 1 - calculate_z(r_outer, p, 1/formal_integral_model.time_explosion), atol=1e-5)
    

@cuda.jit
def populate_z_photosphere_caller(r_inner, r_outer, time_explosion, p, oz, oshell_id, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.populate_z_cuda(r_inner,
                                                     r_outer,
                                                     time_explosion,
                                                     p, 
                                                     oz, 
                                                     oshell_id)


@pytest.mark.parametrize("p", [1e-5, 0.5, 0.99, 1])
def test_populate_z_shells(formal_integral_model, p):
    """
    Test the case where p > r[0]
    
    oz is redshift
    """
    actual = np.zeros(4)
    
    #func = formal_integral_cuda.populate_z_cuda

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
    expected_oshell_id[:expected_N] = (
        np.abs(np.arange(0.5, expected_N, 1) - offset) - 0.5 + idx
    )
    
    expected_oz[0:offset] = 1 + calculate_z(
        r_outer[np.arange(size, idx, -1) - 1], p, inv_t
    )
    expected_oz[offset:expected_N] = 1 - calculate_z(
        r_outer[np.arange(idx, size, 1)], p, inv_t
    )

    
    populate_z_shells_caller[1,4](r_inner, 
                                  r_outer, 
                                  formal_integral_model.time_explosion, 
                                  p, 
                                  oz, 
                                  oshell_id, 
                                  actual)

    #This is to find which p is being tested, so the assert statement will not crash as
    #actual is an array and needs to be indexed
    
    #Need to manually do the math here to see where the error is coming from!
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
    r = 1.0
    func = formal_integral_cuda.calculate_p_values

    expected = r / (N - 1) * np.arange(0, N, dtype=np.float64)
    actual = np.zeros_like(expected, dtype=np.float64)

    actual[::] = func(r, N)
    ntest.assert_allclose(actual, expected)
    ntest.assert_almost_equal(actual, expected)


#Make test for every function call, even the convenience ones

#Review this test with Jack!

@pytest.mark.parametrize(
    ["nu", "nu_insert"],
    [
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 0),
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 1),
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 5),
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 6),
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 7),
        ([0, 1, 1, 2, 3, 4, 5, 6, 6, 7, 8 , 9 , 10], 10),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 4), #This test case will call the else statement
    ],
)
def test_line_search_cuda(nu, nu_insert):
    actual = np.zeros(7)
    nu = np.asarray(nu)
    imin = 0 #Lowest Index
    imax = 12 #Highest index in the input array 
    if nu_insert > nu[imin]:
        expected = imin
    elif nu_insert < nu[imax]:
        expected = imax + 1
    else:
        expected = reverse_binary_search(nu, nu_insert, imin, imax)
        expected = expected + 1
    
    line_search_cuda_caller[1, 7](nu, nu_insert, actual)
    
    nu_loc = 0
    if nu_insert == 0:
        nu_loc = 0
    elif nu_insert == 1:
        nu_loc = 1
    elif nu_insert == 5:
        nu_loc = 2
    elif nu_insert == 6:
        nu_loc = 3
    elif nu_insert == 7:
        nu_loc = 4
    elif nu_insert == 10:
        nu_loc = 5
    elif nu_insert == 4:
        nu_loc = 6
    
    ntest.assert_almost_equal(actual[nu_loc], expected)

def reverse_binary_search(x, x_insert, imin, imax):
    """
    Look for a place to insert a value in an inversely sorted float array.

    Inputs:
        :x: (array) an inversely (largest to lowest) sorted float array
        :x_insert: (value) a value to insert
        :imin: (int) lower bound
        :imax: (int) upper bound

    Outputs:
        index of the next boundary to the left
    """
    # ret_val = TARDIS_ERROR_OK # check
    if (x_insert > x[imin]) or (x_insert < x[imax]):
        raise BoundsError  # check
    return len(x) - 1 - np.searchsorted(x[::-1], x_insert, side="right")

@cuda.jit
def line_search_cuda_caller(nu, nu_insert, actual):
    x = cuda.grid(1)
    actual[x] = formal_integral_cuda.line_search_cuda(nu, nu_insert, 13)
    


@pytest.mark.parametrize(
    ["nu", "nu_insert"],
    [
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 0),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 1),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 5),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 6),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 7),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 10),
        ([10, 9, 8, 7, 6, 6, 5, 4, 3, 2, 1, 1, 0], 4), #This test case will call the else statement
    ],
)
def test_reverse_binary_search(nu, nu_insert):
    actual = np.zeros(7)
    
    nu = np.asarray(nu)
    if (nu_insert > nu[0]) or (nu_insert < nu[12]):
        raise BoundsError  # check
    expected =  len(nu) - 1 - np.searchsorted(nu[::-1], nu_insert, side="right")
    
    actual = np.zeros(7)
    
    cuda_searchsorted_value_right_caller[1, 7](nu, nu_insert, actual)
    
    nu_loc = 0
    if nu_insert == 0:
        nu_loc = 0
    elif nu_insert == 1:
        nu_loc = 1
    elif nu_insert == 5:
        nu_loc = 2
    elif nu_insert == 6:
        nu_loc = 3
    elif nu_insert == 7:
        nu_loc = 4
    elif nu_insert == 10:
        nu_loc = 5
    elif nu_insert == 4:
        nu_loc = 6
    
    ntest.assert_almost_equal(actual[nu_loc], expected)

@cuda.jit
def cuda_searchsorted_value_right_caller(nu, nu_insert, actual):
    x = cuda.grid(1)
    actual[x] = len(nu) - 1 - cuda_searchsorted_value_right(nu[::-1], nu_insert)


