import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry


@pytest.fixture(scope="function")
def homologous_radial1d_geometry():
    velocity = np.arange(8000, 21000, 1000) * u.km / u.s
    v_inner = velocity[:-1]
    v_outer = velocity[1:]
    time_explosion = 5 * u.day
    geometry = HomologousRadial1DGeometry(
        v_inner, v_outer, v_inner[0], v_outer[-1], time_explosion
    )
    return geometry


def test_vb_indices(homologous_radial1d_geometry):
    # Testing if the indices returned are correct when inner and outer
    # boundary are on the innermost and outermost shell

    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[0]
    )
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-1]
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    assert homologous_radial1d_geometry.v_outer_boundary_index == len(
        homologous_radial1d_geometry.v_inner
    )
    vib_index = homologous_radial1d_geometry.v_inner_boundary_index
    vob_index = homologous_radial1d_geometry.v_outer_boundary_index
    assert np.all(
        homologous_radial1d_geometry.v_inner[vib_index:vob_index]
        == homologous_radial1d_geometry.v_inner
    )
    EPSILON_VELOCITY_SHIFT = 1 * u.km / u.s
    # pivoting around the inner boundary of the simulation

    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[0] + EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[0] - EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0

    # pivoting around the first shell boundary of the simulation
    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[1] - EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[1] + EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 1

    # pivoting around the outer boundary of the simulation
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-1] + EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-1] - EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12

    # pivoting around the second to outer boundary of the simulation
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-2] + EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-2] - EPSILON_VELOCITY_SHIFT
    )
    assert homologous_radial1d_geometry.v_outer_boundary_index == 11


def test_velocity_boundary(homologous_radial1d_geometry):
    # testing the active cell boundaries when setting the boundaries

    homologous_radial1d_geometry.v_inner_boundary = 7999 * u.km / u.s
    npt.assert_almost_equal(
        homologous_radial1d_geometry.v_inner_active[0].value, 7999
    )
    assert len(homologous_radial1d_geometry.v_inner_active) == len(
        homologous_radial1d_geometry.v_inner
    )

    homologous_radial1d_geometry.v_inner_boundary = 8001 * u.km / u.s
    npt.assert_almost_equal(
        homologous_radial1d_geometry.v_inner_active[0].value, 8001
    )
    assert len(homologous_radial1d_geometry.v_inner_active) == len(
        homologous_radial1d_geometry.v_inner
    )

    homologous_radial1d_geometry.v_inner_boundary = 9001 * u.km / u.s
    npt.assert_almost_equal(
        homologous_radial1d_geometry.v_inner_active[0].value, 9001
    )
    assert len(homologous_radial1d_geometry.v_inner_active) == (
        len(homologous_radial1d_geometry.v_inner) - 1
    )
    homologous_radial1d_geometry.v_inner_boundary = 9000 * u.km / u.s
    npt.assert_almost_equal(
        homologous_radial1d_geometry.v_inner_active[0].value, 9000
    )
    assert len(homologous_radial1d_geometry.v_inner_active) == (
        len(homologous_radial1d_geometry.v_inner) - 1
    )


def test_v_middle_active_default_boundaries(homologous_radial1d_geometry):
    """Test v_middle_active with default (full range) boundaries"""
    v_middle_active = homologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length
    assert len(v_middle_active) == len(
        homologous_radial1d_geometry.v_inner_active
    )
    assert len(v_middle_active) == len(
        homologous_radial1d_geometry.v_outer_active
    )

    # Check that v_middle_active is the average of v_inner_active and v_outer_active
    expected_v_middle = (
        homologous_radial1d_geometry.v_inner_active
        + homologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # Check correct units
    assert v_middle_active.unit == u.km / u.s


def test_v_middle_active_modified_inner_boundary(homologous_radial1d_geometry):
    """Test v_middle_active when inner boundary is modified"""
    # Set inner boundary to a custom value
    homologous_radial1d_geometry.v_inner_boundary = 9500 * u.km / u.s

    v_middle_active = homologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length (one less shell)
    assert len(v_middle_active) == len(
        homologous_radial1d_geometry.v_inner_active
    )

    # Check computation
    expected_v_middle = (
        homologous_radial1d_geometry.v_inner_active
        + homologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # First middle velocity should be between modified boundary and corresponding outer velocity
    npt.assert_almost_equal(
        v_middle_active[0].value,
        (9500 + 10000) / 2.0,
    )


def test_v_middle_active_modified_outer_boundary(homologous_radial1d_geometry):
    """Test v_middle_active when outer boundary is modified"""
    # Set outer boundary to a custom value
    homologous_radial1d_geometry.v_outer_boundary = 19500 * u.km / u.s

    v_middle_active = homologous_radial1d_geometry.v_middle_active

    # Check that v_middle_active has correct length (one less shell)
    assert len(v_middle_active) == len(
        homologous_radial1d_geometry.v_inner_active
    )

    # Check computation
    expected_v_middle = (
        homologous_radial1d_geometry.v_inner_active
        + homologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )

    # Last middle velocity should be between corresponding inner velocity and modified boundary
    npt.assert_almost_equal(
        v_middle_active[-1].value,
        (19000 + 19500) / 2.0,
    )


def test_v_middle_active_both_boundaries_modified(homologous_radial1d_geometry):
    """Test v_middle_active when both boundaries are modified"""
    # Set both boundaries to custom values
    homologous_radial1d_geometry.v_inner_boundary = 9500 * u.km / u.s
    homologous_radial1d_geometry.v_outer_boundary = 19500 * u.km / u.s

    v_middle_active = homologous_radial1d_geometry.v_middle_active

    # Check computation
    expected_v_middle = (
        homologous_radial1d_geometry.v_inner_active
        + homologous_radial1d_geometry.v_outer_active
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


def test_v_middle_active_narrow_range(homologous_radial1d_geometry):
    """Test v_middle_active with a narrow range of active shells"""
    # Set boundaries to include only a few shells
    homologous_radial1d_geometry.v_inner_boundary = 11500 * u.km / u.s
    homologous_radial1d_geometry.v_outer_boundary = 14500 * u.km / u.s

    v_middle_active = homologous_radial1d_geometry.v_middle_active

    # Should have 4 active shells: [11500-12000] and [12000-13000] and [13000-14000] and [14000-14500]
    assert len(v_middle_active) == 4

    # Check computation
    expected_v_middle = (
        homologous_radial1d_geometry.v_inner_active
        + homologous_radial1d_geometry.v_outer_active
    ) / 2.0
    npt.assert_array_almost_equal(
        v_middle_active.value, expected_v_middle.value
    )
