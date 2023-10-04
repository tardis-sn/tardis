from astropy import units as u
import numpy as np

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


def test_v_indices(homologous_radial1d_geometry):
    # Testing if the indices returned are correct when inner and outer
    # boundary are on the innermost and outermost shell

    homologous_radial1d_geometry.v_inner_boundary = (
        homologous_radial1d_geometry.v_inner[0]
    )
    homologous_radial1d_geometry.v_outer_boundary = (
        homologous_radial1d_geometry.v_outer[-1]
    )
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12
    vib_index = homologous_radial1d_geometry.v_inner_boundary_index
    vob_index = homologous_radial1d_geometry.v_outer_boundary_index
    assert np.all(
        homologous_radial1d_geometry.v_inner[vib_index:vob_index]
        == homologous_radial1d_geometry.v_inner
    )

    #pivoting around the inner boundary of the simulation

    homologous_radial1d_geometry.v_inner_boundary = homologous_radial1d_geometry.v_inner[0] + 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    homologous_radial1d_geometry.v_inner_boundary = homologous_radial1d_geometry.v_inner[0] - 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0

    #pivoting around the first shell boundary of the simulation
    homologous_radial1d_geometry.v_inner_boundary = homologous_radial1d_geometry.v_inner[1] - 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_inner_boundary_index == 0
    homologous_radial1d_geometry.v_inner_boundary = homologous_radial1d_geometry.v_inner[1] + 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_inner_boundary_index == 1

    #pivoting around the outer boundary of the simulation
    homologous_radial1d_geometry.v_outer_boundary = homologous_radial1d_geometry.v_outer[-1] + 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12
    homologous_radial1d_geometry.v_outer_boundary = homologous_radial1d_geometry.v_outer[-1] - 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12

    #pivoting around the second to outer boundary of the simulation
    homologous_radial1d_geometry.v_outer_boundary = homologous_radial1d_geometry.v_outer[-2] + 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_outer_boundary_index == 12
    homologous_radial1d_geometry.v_outer_boundary = homologous_radial1d_geometry.v_outer[-2] - 0.0001 * u.km / u.s
    assert homologous_radial1d_geometry.v_outer_boundary_index == 11

