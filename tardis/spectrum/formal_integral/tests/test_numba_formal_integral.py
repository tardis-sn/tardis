import numpy as np
import numpy.testing as ntest
import pytest

from tardis import constants as c
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.spectrum.formal_integral.base import C_INV
from tardis.spectrum.formal_integral.formal_integral import FormalIntegrator
import tardis.spectrum.formal_integral.formal_integral_numba as formal_integral_numba


TESTDATA = [
    {
        "r": np.linspace(1, 2, 3, dtype=np.float64),
    },
    {
        "r": np.linspace(0, 1, 3),
    },
    # {"r": np.linspace(1, 2, 10, dtype=np.float64)},
]
@pytest.fixture(scope="function", params=TESTDATA)
def formal_integral_geometry(request):
    r = request.param["r"]
    geometry = NumbaRadial1DGeometry(
        r[:-1],
        r[1:],
        r[:-1] * c.c.cgs.value,
        r[1:] * c.c.cgs.value,
    )
    return geometry


@pytest.fixture(scope="function")
def time_explosion():
    return 1 / c.c.cgs.value


@pytest.mark.parametrize("N", (1e2, 1e3, 1e4, 1e5))
def test_trapezoid_integration(N):
    h = 1.0
    N = int(N)
    data = np.random.random(N)

    actual = formal_integral_numba.trapezoid_integration(data, h)
    expected = np.trapz(data)

    ntest.assert_almost_equal(actual, expected)


def calculate_z(r, p):
    return np.sqrt(r * r - p * p)


@pytest.mark.parametrize("p", [0.0, 0.5, 1.0])
def test_calculate_z(formal_integral_geometry, time_explosion, p):

    inv_t = 1.0 / time_explosion
    size = len(formal_integral_geometry.r_outer)
    r_outer = formal_integral_geometry.r_outer
    for r in r_outer:

        actual = formal_integral_numba.calculate_z(r, p, inv_t)
        if p >= r:
            assert actual == 0
        else:
            desired = np.sqrt(r * r - p * p) * C_INV * inv_t
            ntest.assert_almost_equal(actual, desired)


@pytest.mark.parametrize("p", [0, 0.5, 1])
def test_populate_z_photosphere(formal_integral_geometry, time_explosion, p):
    """
    Test the case where p < r[0]
    That means we 'hit' all shells from inside to outside.
    """
    size = len(formal_integral_geometry.r_outer)
    r_inner = formal_integral_geometry.r_inner
    r_outer = formal_integral_geometry.r_outer

    p = r_inner[0] * p
    oz = np.zeros_like(r_inner)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    N = formal_integral_numba.populate_z(formal_integral_geometry, time_explosion, p, oz, oshell_id)
    assert N == size

    ntest.assert_allclose(oshell_id, np.arange(0, size, 1))

    ntest.assert_allclose(oz, 1 - calculate_z(r_outer, p), atol=1e-5)


@pytest.mark.parametrize("p", [1e-5, 0.5, 0.99, 1])
def test_populate_z_shells(formal_integral_geometry, time_explosion, p):
    """
    Test the case where p > r[0]
    """

    size = len(formal_integral_geometry.r_inner)
    r_inner = formal_integral_geometry.r_inner
    r_outer = formal_integral_geometry.r_outer

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

    N = formal_integral_numba.populate_z(formal_integral_geometry, time_explosion, p, oz, oshell_id)

    assert N == expected_N

    ntest.assert_allclose(oshell_id, expected_oshell_id)

    ntest.assert_allclose(oz, expected_oz, atol=1e-5)
