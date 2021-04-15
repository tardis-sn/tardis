import pytest
import numpy.testing as npt
import numpy as np

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

    x_old = gamma_ray.location.x
    y_old = gamma_ray.location.y
    z_old = gamma_ray.location.z

    x_new = x_old + distance * gamma_ray.direction.x
    y_new = y_old + distance * gamma_ray.direction.y
    z_new = z_old + distance * gamma_ray.direction.z

    expected_r = np.sqrt(x_new ** 2.0 + y_new ** 2.0 + z_new ** 2.0)
    expected_mu = z_new / expected_r

    actual = move_gamma_ray(gamma_ray, distance)

    npt.assert_almost_equal(actual.location.r, expected_r)
    npt.assert_almost_equal(actual.location.mu, expected_mu)


@pytest.mark.xfail(reason="To be implemented")
def test_density_sampler():
    assert False
