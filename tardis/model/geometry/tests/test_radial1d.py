from astropy import units as u
import numpy as np
import numpy.testing as npt

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry

import pytest


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
