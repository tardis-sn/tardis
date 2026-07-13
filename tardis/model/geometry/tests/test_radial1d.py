import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.model.geometry.radial1d_nonhomologous import (
    NonhomologousRadial1DGeometry,
    NumbaNonhomologousRadial1DGeometry,
)


@pytest.fixture(scope="function")
def nonhomologous_radial1d_geometry():
    velocity = np.arange(8000, 21000, 1000) * u.km / u.s
    v_inner = velocity[:-1]
    v_outer = velocity[1:]
    time_explosion = 5 * u.day
    return NonhomologousRadial1DGeometry(
        r_inner=(v_inner * time_explosion).cgs,
        r_outer=(v_outer * time_explosion).cgs,
        v_inner=v_inner,
        v_outer=v_outer,
        r_inner_boundary=(v_inner[0] * time_explosion).cgs,
        r_outer_boundary=(v_outer[-1] * time_explosion).cgs,
        v_inner_boundary=v_inner[0],
        v_outer_boundary=v_outer[-1],
    )


def test_vb_indices(nonhomologous_radial1d_geometry):
    # Testing if the indices returned are correct when inner and outer
    # boundary are on the innermost and outermost shell

    nonhomologous_radial1d_geometry.v_inner_boundary = (
        nonhomologous_radial1d_geometry.v_inner[0]
    )
    nonhomologous_radial1d_geometry.v_outer_boundary = (
        nonhomologous_radial1d_geometry.v_outer[-1]
    )
    assert nonhomologous_radial1d_geometry.v_inner_boundary_index == 0
    assert nonhomologous_radial1d_geometry.v_outer_boundary_index == len(
        nonhomologous_radial1d_geometry.v_inner
    )
    vib_index = nonhomologous_radial1d_geometry.v_inner_boundary_index
    vob_index = nonhomologous_radial1d_geometry.v_outer_boundary_index
    assert np.all(
        nonhomologous_radial1d_geometry.v_inner[vib_index:vob_index]
        == nonhomologous_radial1d_geometry.v_inner
    )
    EPSILON_VELOCITY_SHIFT = 1 * u.km / u.s
    # pivoting around the inner boundary of the simulation

    nonhomologous_radial1d_geometry.v_inner_boundary = (
        nonhomologous_radial1d_geometry.v_inner[0] + EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_inner_boundary_index == 0
    nonhomologous_radial1d_geometry.v_inner_boundary = (
        nonhomologous_radial1d_geometry.v_inner[0] - EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_inner_boundary_index == 0

    # pivoting around the first shell boundary of the simulation
    nonhomologous_radial1d_geometry.v_inner_boundary = (
        nonhomologous_radial1d_geometry.v_inner[1] - EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_inner_boundary_index == 0
    nonhomologous_radial1d_geometry.v_inner_boundary = (
        nonhomologous_radial1d_geometry.v_inner[1] + EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_inner_boundary_index == 1

    # pivoting around the outer boundary of the simulation
    nonhomologous_radial1d_geometry.v_outer_boundary = (
        nonhomologous_radial1d_geometry.v_outer[-1] + EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_outer_boundary_index == 12
    nonhomologous_radial1d_geometry.v_outer_boundary = (
        nonhomologous_radial1d_geometry.v_outer[-1] - EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_outer_boundary_index == 12

    # pivoting around the second to outer boundary of the simulation
    nonhomologous_radial1d_geometry.v_outer_boundary = (
        nonhomologous_radial1d_geometry.v_outer[-2] + EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_outer_boundary_index == 12
    nonhomologous_radial1d_geometry.v_outer_boundary = (
        nonhomologous_radial1d_geometry.v_outer[-2] - EPSILON_VELOCITY_SHIFT
    )
    assert nonhomologous_radial1d_geometry.v_outer_boundary_index == 11


def test_numba_nonhomologous_velocity(nonhomologous_radial1d_geometry):
    numba_geometry = nonhomologous_radial1d_geometry.to_numba()
    radius = numba_geometry.r_inner[0] * 1.2

    npt.assert_allclose(
        numba_geometry.velocity_gradient[0],
        (numba_geometry.v_outer[0] - numba_geometry.v_inner[0])
        / (numba_geometry.r_outer[0] - numba_geometry.r_inner[0]),
    )
    npt.assert_allclose(
        numba_geometry.get_velocity(radius, 0),
        numba_geometry.v_inner[0]
        + numba_geometry.velocity_gradient[0]
        * (radius - numba_geometry.r_inner[0]),
    )


def test_numba_nonhomologous_velocity_at_origin():
    time_explosion = 5.0
    r_inner = np.array([0.0, 10.0])
    r_outer = np.array([10.0, 20.0])
    v_inner = r_inner / time_explosion
    v_outer = r_outer / time_explosion

    numba_geometry = NumbaNonhomologousRadial1DGeometry(
        r_inner,
        r_outer,
        v_inner,
        v_outer,
    )

    npt.assert_allclose(numba_geometry.get_velocity(15.0, 1), 3.0)


def test_velocity_boundary(nonhomologous_radial1d_geometry):
    # testing the active cell boundaries when setting the boundaries

    nonhomologous_radial1d_geometry.v_inner_boundary = 7999 * u.km / u.s
    npt.assert_almost_equal(
        nonhomologous_radial1d_geometry.v_inner_active[0].value, 7999
    )
    assert len(nonhomologous_radial1d_geometry.v_inner_active) == len(
        nonhomologous_radial1d_geometry.v_inner
    )

    nonhomologous_radial1d_geometry.v_inner_boundary = 8001 * u.km / u.s
    npt.assert_almost_equal(
        nonhomologous_radial1d_geometry.v_inner_active[0].value, 8001
    )
    assert len(nonhomologous_radial1d_geometry.v_inner_active) == len(
        nonhomologous_radial1d_geometry.v_inner
    )

    nonhomologous_radial1d_geometry.v_inner_boundary = 9001 * u.km / u.s
    npt.assert_almost_equal(
        nonhomologous_radial1d_geometry.v_inner_active[0].value, 9001
    )
    assert len(nonhomologous_radial1d_geometry.v_inner_active) == (
        len(nonhomologous_radial1d_geometry.v_inner) - 1
    )
    nonhomologous_radial1d_geometry.v_inner_boundary = 9000 * u.km / u.s
    npt.assert_almost_equal(
        nonhomologous_radial1d_geometry.v_inner_active[0].value, 9000
    )
    assert len(nonhomologous_radial1d_geometry.v_inner_active) == (
        len(nonhomologous_radial1d_geometry.v_inner) - 1
    )


def test_velocity_boundary_updates_radius_boundary(
    nonhomologous_radial1d_geometry,
):
    time_explosion = 5 * u.day
    v_inner_boundary = 9500 * u.km / u.s
    v_outer_boundary = 19500 * u.km / u.s

    nonhomologous_radial1d_geometry.v_inner_boundary = v_inner_boundary
    nonhomologous_radial1d_geometry.v_outer_boundary = v_outer_boundary

    assert u.isclose(
        nonhomologous_radial1d_geometry.r_inner_active[0],
        (v_inner_boundary * time_explosion).cgs,
    )
    assert u.isclose(
        nonhomologous_radial1d_geometry.r_outer_active[-1],
        (v_outer_boundary * time_explosion).cgs,
    )


def test_v_middle_active_default_boundaries(nonhomologous_radial1d_geometry):
    """Test v_middle_active with default (full range) boundaries"""
    v_middle_active = nonhomologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length
    assert len(v_middle_active) == len(
        nonhomologous_radial1d_geometry.v_inner_active
    )
    assert len(v_middle_active) == len(
        nonhomologous_radial1d_geometry.v_outer_active
    )

    # Check that v_middle_active is the average of v_inner_active and v_outer_active
    expected_v_middle = (
        nonhomologous_radial1d_geometry.v_inner_active
        + nonhomologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # Check correct units
    assert v_middle_active.unit == u.km / u.s


def test_v_middle_active_modified_inner_boundary(nonhomologous_radial1d_geometry):
    """Test v_middle_active when inner boundary is modified"""
    # Set inner boundary to a custom value
    nonhomologous_radial1d_geometry.v_inner_boundary = 9500 * u.km / u.s

    v_middle_active = nonhomologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length (one less shell)
    assert len(v_middle_active) == len(
        nonhomologous_radial1d_geometry.v_inner_active
    )

    # Check computation
    expected_v_middle = (
        nonhomologous_radial1d_geometry.v_inner_active
        + nonhomologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # First middle velocity should be between modified boundary and corresponding outer velocity
    npt.assert_almost_equal(
        v_middle_active[0].value,
        (9500 + 10000) / 2.0,
    )


def test_v_middle_active_modified_outer_boundary(nonhomologous_radial1d_geometry):
    """Test v_middle_active when outer boundary is modified"""
    # Set outer boundary to a custom value
    nonhomologous_radial1d_geometry.v_outer_boundary = 19500 * u.km / u.s

    v_middle_active = nonhomologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length (one less shell)
    assert len(v_middle_active) == len(
        nonhomologous_radial1d_geometry.v_inner_active
    )

    # Check computation
    expected_v_middle = (
        nonhomologous_radial1d_geometry.v_inner_active
        + nonhomologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # Last middle velocity should be between corresponding inner velocity and modified boundary
    npt.assert_almost_equal(
        v_middle_active[-1].value,
        (19000 + 19500) / 2.0,
    )


def test_v_middle_active_both_boundaries_modified(nonhomologous_radial1d_geometry):
    """Test v_middle_active when both boundaries are modified"""
    # Set both boundaries to custom values
    nonhomologous_radial1d_geometry.v_inner_boundary = 9500 * u.km / u.s
    nonhomologous_radial1d_geometry.v_outer_boundary = 19500 * u.km / u.s

    v_middle_active = nonhomologous_radial1d_geometry.v_middle_active

    # Check computation
    expected_v_middle = (
        nonhomologous_radial1d_geometry.v_inner_active
        + nonhomologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # Check first and last values
    npt.assert_almost_equal(
        v_middle_active[0].value,
        (9500 + 10000) / 2.0,
    )
    npt.assert_almost_equal(
        v_middle_active[-1].value,
        (19000 + 19500) / 2.0,
    )


def test_v_middle_active_narrow_range(nonhomologous_radial1d_geometry):
    """Test v_middle_active with a narrow range of active shells"""
    # Set boundaries to include only a few shells
    nonhomologous_radial1d_geometry.v_inner_boundary = 11500 * u.km / u.s
    nonhomologous_radial1d_geometry.v_outer_boundary = 14500 * u.km / u.s

    v_middle_active = nonhomologous_radial1d_geometry.v_middle_active

    # Should have 4 active shells: [11500-12000] and [12000-13000] and [13000-14000] and [14000-14500]
    assert len(v_middle_active) == 4

    # Check computation
    expected_v_middle = (
        nonhomologous_radial1d_geometry.v_inner_active
        + nonhomologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )
