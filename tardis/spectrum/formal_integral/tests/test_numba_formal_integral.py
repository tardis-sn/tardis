import numpy as np
import numpy.testing as ntest
import pytest
from astropy import units as u

import tardis.spectrum.formal_integral.formal_integral_numba as formal_integral_numba
from tardis import constants as c
from tardis.model.geometry.radial1d_homologous import HomologousRadial1DGeometry
from tardis.spectrum.formal_integral.base import C_INV

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
    time_explosion = (1 / c.c.cgs.value) * u.s
    return HomologousRadial1DGeometry(
        r[:-1] * u.cm / time_explosion,
        r[1:] * u.cm / time_explosion,
        None,
        None,
        time_explosion,
    ).to_numba()


@pytest.fixture(scope="function")
def time_explosion():
    return 1 / c.c.cgs.value


def calculate_intersection_point(r, p):
    return np.sqrt(r * r - p * p)


@pytest.mark.parametrize("p", [0.0, 0.5, 1.0])
def test_calculate_intersection_point(formal_integral_geometry, time_explosion, p):
    inv_t = 1.0 / time_explosion
    size = len(formal_integral_geometry.r_outer)
    r_outer = formal_integral_geometry.r_outer
    for r in r_outer:
        actual = formal_integral_numba.calculate_intersection_point(r, p, inv_t)
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
    oz = np.zeros(size + 1)
    oshell_id = np.zeros_like(oz, dtype=np.int64)

    N = formal_integral_numba.populate_intersection_points(
        formal_integral_geometry, time_explosion, p, oz, oshell_id
    )
    assert N == size + 1

    expected_oshell_id = np.minimum(np.arange(size + 1), size - 1)
    ntest.assert_allclose(oshell_id, expected_oshell_id)

    expected_oz = 1 - calculate_intersection_point(
        np.concatenate([r_inner[:1], r_outer]), p
    )
    ntest.assert_allclose(oz, expected_oz, atol=1e-5)


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

    # Each ID is the shell traversed by the segment starting at that point:
    # inwards through shells size - 1 ... idx on the far side, then outwards
    # through shells idx + 1 ... size - 1 on the near side.
    expected_oshell_id[:offset] = np.arange(size - 1, idx - 1, -1)
    expected_oshell_id[offset:expected_N] = np.minimum(
        np.arange(idx + 1, size + 1), size - 1
    )

    expected_oz[0:offset] = 1 + calculate_intersection_point(
        r_outer[np.arange(size, idx, -1) - 1], p
    )
    expected_oz[offset:expected_N] = 1 - calculate_intersection_point(
        r_outer[np.arange(idx, size, 1)], p
    )

    N = formal_integral_numba.populate_intersection_points(
        formal_integral_geometry, time_explosion, p, oz, oshell_id
    )

    assert N == expected_N

    ntest.assert_allclose(oshell_id, expected_oshell_id)

    ntest.assert_allclose(oz, expected_oz, atol=1e-5)


@pytest.mark.parametrize("p", [0.0, 0.5, 0.99, 1.0, 1.05, 1.3, 1.8])
def test_populate_intersection_points_segment_shells(time_explosion, p):
    """
    Test that every ray segment is assigned the shell it actually lies in,
    on both the far and the near side of the ray.
    """
    r = np.linspace(1, 2, 6)
    geometry = HomologousRadial1DGeometry(
        r[:-1] * u.cm / (time_explosion * u.s),
        r[1:] * u.cm / (time_explosion * u.s),
        None,
        None,
        time_explosion * u.s,
    ).to_numba()
    size = len(geometry.r_inner)

    oz = np.zeros(2 * size)
    oshell_id = np.zeros_like(oz, dtype=np.int64)
    N = formal_integral_numba.populate_intersection_points(
        geometry, time_explosion, p, oz, oshell_id
    )

    # with c * time_explosion = 1, the line-of-sight coordinate is 1 - z
    segment_midpoints = 1 - 0.5 * (oz[: N - 1] + oz[1:N])
    segment_radii = np.sqrt(p**2 + segment_midpoints**2)
    segment_shells = oshell_id[: N - 1]

    assert np.all(segment_radii >= geometry.r_inner[segment_shells])
    assert np.all(segment_radii <= geometry.r_outer[segment_shells])
    if p <= geometry.r_inner[0]:
        # rays hitting the photosphere start at the photosphere
        ntest.assert_allclose(
            oz[0],
            1 - calculate_intersection_point(geometry.r_inner[0], p),
            atol=1e-5,
        )
