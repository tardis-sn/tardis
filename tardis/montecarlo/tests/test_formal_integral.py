import pytest
import numpy as np
from tardis import constants as c

from numpy.ctypeslib import (
    as_array,
    as_ctypes,
)


import numpy.testing as ntest

from tardis.util.base import intensity_black_body
import tardis.montecarlo.formal_integral as formal_integral


@pytest.mark.parametrize(
    ["nu", "T"],
    [
        (1e14, 1e4),
        (0, 1),
        (1, 1),
    ],
)
def test_intensity_black_body(nu, T):
    func = formal_integral.intensity_black_body
    actual = func(nu, T)
    print(actual, type(actual))
    expected = intensity_black_body(nu, T)
    ntest.assert_almost_equal(actual, expected)


@pytest.mark.parametrize("N", (1e2, 1e3, 1e4, 1e5))
def test_trapezoid_integration(N):
    func = formal_integral.trapezoid_integration
    h = 1.0
    N = int(N)
    data = np.random.random(N)

    actual = func(data, h, int(N))
    expected = np.trapz(data)

    ntest.assert_almost_equal(actual, expected)


@pytest.mark.skipif(
    True, reason="static inline functions are not inside the library"
)
def test_calculate_z():
    pass


def calculate_z(r, p):
    return np.sqrt(r * r - p * p)


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
def formal_integral_model(request, model):
    r = request.param["r"]
    model.no_of_shells_i = r.shape[0] - 1
    model.time_explosion = 1 / c.c.cgs.value
    model.r_outer_i.contents = as_ctypes(r[1:])
    model.r_inner_i.contents = as_ctypes(r[:-1])
    return model


@pytest.mark.xfail(reason="not implemented")
@pytest.mark.parametrize("p", [0, 0.5, 1])
def test_populate_z_photosphere(formal_integral_model, p):
    """
    Test the case where p < r[0]
    That means we 'hit' all shells from inside to outside.
    """
    integrator = formal_integral.FormalIntegrator(
        formal_integral_model, None, None
    )
    func = integrator.populate_z
    size = formal_integral_model.no_of_shells_i
    r_inner = as_array(formal_integral_model.r_inner_i, (size,))
    r_outer = as_array(formal_integral_model.r_outer_i, (size,))

    p = r_inner[0] * p
    oz = np.zeros_like(r_inner)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    N = func(p, oz, oshell_id)
    assert N == size

    ntest.assert_allclose(oshell_id, np.arange(0, size, 1))

    ntest.assert_allclose(oz, 1 - calculate_z(r_outer, p), atol=1e-5)


@pytest.mark.parametrize("p", [1e-5, 0.5, 0.99, 1])
def test_populate_z_shells(formal_integral_model, p):
    """
    Test the case where p > r[0]
    """
    integrator = formal_integral.FormalIntegrator(
        formal_integral_model, None, None
    )
    func = integrator.populate_z

    size = formal_integral_model.no_of_shells_i
    r_inner = as_array(formal_integral_model.r_inner_i, (size,))
    r_outer = as_array(formal_integral_model.r_outer_i, (size,))

    p = r_inner[0] + (r_outer[-1] - r_inner[0]) * p
    idx = np.searchsorted(r_outer, p, side="right")

    oz = np.zeros(size * 2)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    offset = size - idx

    expected_N = (offset) * 2
    expected_oz = np.zeros_like(oz)
    expected_oshell_id = np.zeros_like(oshell_id)

    # Calculated way to determine which shells get hit
    expected_oshell_id[:expected_N] = (
        np.abs(np.arange(0.5, expected_N, 1) - offset) - 0.5 + idx
    )

    expected_oz[0:offset] = 1 + calculate_z(
        r_outer[np.arange(size, idx, -1) - 1], p
    )
    expected_oz[offset:expected_N] = 1 - calculate_z(
        r_outer[np.arange(idx, size, 1)], p
    )

    N = func(p, oz, oshell_id)

    assert N == expected_N

    ntest.assert_allclose(oshell_id, expected_oshell_id)

    ntest.assert_allclose(oz, expected_oz, atol=1e-5)


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
    func = formal_integral.calculate_p_values

    expected = r / (N - 1) * np.arange(0, N, dtype=np.float64)
    actual = np.zeros_like(expected, dtype=np.float64)

    func(r, N, actual)
    ntest.assert_allclose(actual, expected)
