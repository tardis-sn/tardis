import pytest
import numpy.testing as npt
import numpy as np
from astropy.coordinates import cartesian_to_spherical

from tardis.energy_input.gamma_ray_grid import (
    calculate_distance_radial,
    distance_trace,
    move_gamma_ray,
    density_sampler,
)


@pytest.mark.xfail(reason="To be implemented")
def test_calculate_distance_radial():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_distance_trace():
    assert False


def test_move_gamma_ray(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    gamma_ray = basic_gamma_ray
    distance = 1.0e15

    x_old, y_old, z_old = gamma_ray.location.get_cartesian_coords
    x_dir, y_dir, z_dir = gamma_ray.direction.get_cartesian_coords

    x_new = x_old + distance * (1 + 1e-7) * x_dir
    y_new = y_old + distance * (1 + 1e-7) * y_dir
    z_new = z_old + distance * (1 + 1e-7) * z_dir

    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)

    expected_r = r.value
    expected_theta = theta.value + 0.5 * np.pi
    expected_phi = phi.value

    actual = move_gamma_ray(gamma_ray, distance)

    npt.assert_almost_equal(actual.location.r, expected_r)
    npt.assert_almost_equal(actual.location.theta, expected_theta)
    npt.assert_almost_equal(actual.location.phi, expected_phi)


@pytest.mark.xfail(reason="To be implemented")
def test_density_sampler():
    assert False
