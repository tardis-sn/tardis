import pytest
import numpy as np
from tardis import constants as c
from astropy import units as u

from copy import deepcopy
import numpy.testing as ntest
from numba import cuda
from numba import njit
from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel


import tardis.montecarlo.montecarlo_numba.formal_integral_cuda as formal_integral_cuda
import tardis.montecarlo.montecarlo_numba.formal_integral as formal_integral_numba
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.montecarlo.montecarlo_numba.numba_interface import NumbaModel


from tardis.montecarlo.montecarlo_numba.formal_integral import FormalIntegrator
from tardis.montecarlo.montecarlo_numba.formal_integral import (
    NumbaFormalIntegrator,
)

from tardis.montecarlo import MontecarloRunner


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
def test_intensity_black_body_cuda(nu, T):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    intensity_black_body to 15 decimals. This
    is done as both results have 15 digits of precision.
    """
    actual = np.zeros(3)
    black_body_caller[1, 3](nu, T, actual)

    expected = formal_integral_numba.intensity_black_body(nu, T)

    ntest.assert_allclose(actual, expected, rtol=1e-14)


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
@pytest.mark.parametrize("N", (1e2, 1e3, 1e4, 1e5))
def test_trapezoid_integration_cuda(N):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    trapezoid_integration to 15 decimals. This
    is done as both results have 15 digits of precision.
    """
    actual = np.zeros(1)

    h = 1.0
    N = int(N)
    np.random.seed(12)
    data = np.random.random(N)

    expected = formal_integral_numba.trapezoid_integration(data, h)
    trapezoid_integration_caller[1, 1](data, h, actual)

    # This is 1e-13, as more points are added to the integration
    # there will be more floating point error due to the difference
    # in how the trapezoid integration is called.
    ntest.assert_allclose(actual[0], expected, rtol=1e-13)


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
def formal_integral_geometry(request):
    """
    This gets the Numba model to be used in later tests
    """
    r = request.param["r"]
    geometry = NumbaRadial1DGeometry(
        r[:-1],
        r[1:],
        r[:-1] * c.c.cgs.value,
        r[1:] * c.c.cgs.value,
    )
    return geometry


@pytest.fixture(scope="function")
def formal_integral_model():
    """
    This gets the Numba model to be used in later tests
    """
    model = NumbaModel(
        1 / c.c.cgs.value,
    )
    return model


@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize(["p", "p_loc"], [(0.0, 0), (0.5, 1), (1.0, 2)])
def test_calculate_z_cuda(
    formal_integral_geometry, formal_integral_model, p, p_loc
):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    calculate_z to 15 decimals. This is done as
    both results have 15 digits of precision.
    """
    actual = np.zeros(3)
    inv_t = 1.0 / formal_integral_model.time_explosion
    size = len(formal_integral_geometry.r_outer)
    r_outer = formal_integral_geometry.r_outer
    for r in r_outer:
        calculate_z_caller[1, 3](r, p, inv_t, actual)
        expected = formal_integral_numba.calculate_z(r, p, inv_t)

        ntest.assert_allclose(actual[p_loc], expected, rtol=1e-14)


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
def test_populate_z(formal_integral_geometry, formal_integral_model, p, p_loc):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    populate_z
    """
    size = len(formal_integral_geometry.r_inner)
    oz = np.zeros(size * 2)
    expected_oz = np.zeros(size * 2)
    oshell_id = np.zeros_like(oz, dtype=np.int64)
    expected_oshell_id = np.zeros_like(oz, dtype=np.int64)

    expected = formal_integral_numba.populate_z(
        formal_integral_geometry,
        formal_integral_model,
        p,
        expected_oz,
        expected_oshell_id,
    )

    actual = np.zeros(6)
    populate_z_caller[1, 6](
        formal_integral_geometry.r_inner,
        formal_integral_geometry.r_outer,
        formal_integral_model.time_explosion,
        p,
        oz,
        oshell_id,
        actual,
    )

    ntest.assert_equal(actual[p_loc], expected)
    ntest.assert_equal(oshell_id, expected_oshell_id)
    ntest.assert_allclose(oz, expected_oz, atol=1e-4)


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
def test_calculate_p_values(N):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    calculate_p_values to 15 decimals. This is done as
    both results have 15 digits of precision.
    """
    r = 1.0

    expected = formal_integral_numba.calculate_p_values(r, N)

    actual = np.zeros_like(expected, dtype=np.float64)
    actual[::] = formal_integral_cuda.calculate_p_values(r, N)

    ntest.assert_allclose(actual, expected, rtol=1e-14)


@pytest.mark.skipif(
    not GPUs_available, reason="No GPU is available to test CUDA function"
)
@pytest.mark.parametrize("nu_insert", np.linspace(3e12, 3e16, 10))
def test_line_search_cuda(nu_insert, verysimple_numba_plasma):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    line_search
    """
    actual = np.zeros(1)
    expected = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu

    expected[0] = formal_integral_numba.line_search(
        line_list_nu, nu_insert, len(line_list_nu)
    )

    line_search_cuda_caller[1, 1](line_list_nu, nu_insert, actual)

    ntest.assert_equal(actual, expected)


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
def test_reverse_binary_search(nu_insert, verysimple_numba_plasma):
    """
    Initializes the test of the cuda version
    against the numba implementation of the
    reverse_binary_search. The one extra input not included
    in np.linspace a low edge case for testing.
    """
    actual = np.zeros(1)
    expected = np.zeros(1)
    line_list_nu = verysimple_numba_plasma.line_list_nu

    imin = 0
    imax = len(line_list_nu) - 1

    expected[0] = formal_integral_numba.reverse_binary_search(
        line_list_nu, nu_insert, imin, imax
    )
    reverse_binary_search_cuda_caller[1, 1](
        line_list_nu, nu_insert, imin, imax, actual
    )

    ntest.assert_equal(actual, expected)


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
def test_full_formal_integral(
    no_of_packets, iterations, config_verysimple, simulation_verysimple
):
    """
    This function initializes both the cuda and numba formal_integrator,
    and the runs them and compares results to the 15th decimal place.
    This is done as both results have 15 digits of precision.
    """
    sim = simulation_verysimple

    formal_integrator_numba = FormalIntegrator(
        sim.model, sim.plasma, sim.runner
    )

    formal_integrator_cuda = FormalIntegrator(sim.model, sim.plasma, sim.runner)

    # The function calculate_spectrum sets this property, but in order to test the CUDA.
    # version it is done manually, as well as to speed up the test.
    formal_integrator_numba.interpolate_shells = max(
        2 * formal_integrator_numba.model.no_of_shells, 80
    )

    formal_integrator_cuda.interpolate_shells = max(
        2 * formal_integrator_cuda.model.no_of_shells, 80
    )

    res_numba = formal_integrator_numba.make_source_function()
    att_S_ul_numba = res_numba[0].flatten(order="F")
    Jred_lu_numba = res_numba[1].values.flatten(order="F")
    Jblue_lu_numba = res_numba[2].flatten(order="F")

    res_cuda = formal_integrator_cuda.make_source_function()
    att_S_ul_cuda = res_cuda[0].flatten(order="F")
    Jred_lu_cuda = res_cuda[1].values.flatten(order="F")
    Jblue_lu_cuda = res_cuda[2].flatten(order="F")

    formal_integrator_numba.generate_numba_objects()
    # This is to force the for formal_integrator_numba to use the numba version
    # as it is automatically set to the CUDA version when there is a GPU available
    formal_integrator_numba.integrator = NumbaFormalIntegrator(
        formal_integrator_numba.numba_radial_1d_geometry,
        formal_integrator_numba.numba_model,
        formal_integrator_numba.numba_plasma,
        formal_integrator_numba.points,
    )

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
    )[0]

    L_numba = formal_integrator_numba.integrator.formal_integral(
        formal_integrator_numba.model.t_inner,
        sim.runner.spectrum.frequency,
        sim.runner.spectrum.frequency.shape[0],
        att_S_ul_numba,
        Jred_lu_numba,
        Jblue_lu_numba,
        formal_integrator_numba.runner.tau_sobolevs_integ,
        formal_integrator_numba.runner.electron_densities_integ,
        formal_integrator_numba.points,
    )[0]

    ntest.assert_allclose(L_cuda, L_numba, rtol=1e-14)
