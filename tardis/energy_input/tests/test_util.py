import pytest
import astropy.units as u
import numpy.testing as npt
import numpy as np

import tardis.energy_input.util as util
from tardis.energy_input.util import R_ELECTRON_SQUARED
from tardis import constants as const


@pytest.mark.parametrize(
    ["r", "theta", "phi", "expected_x", "expected_y", "expected_z"],
    [
        (1, 0, 0, 0, 0, 1),
        (1, np.pi, 0, 0, 0, -1),
        (1, np.pi / 2, 0, 1, 0, 0),
        (1, np.pi / 2, np.pi, -1, 0, 0),
        (1, np.pi / 2, np.pi / 2, 0, 1, 0),
    ],
)
def test_spherical_to_cartesian(
    r, theta, phi, expected_x, expected_y, expected_z
):
    actual_x, actual_y, actual_z = util.spherical_to_cartesian(r, theta, phi)
    npt.assert_almost_equal(actual_x, expected_x)
    npt.assert_almost_equal(actual_y, expected_y)
    npt.assert_almost_equal(actual_z, expected_z)


@pytest.mark.xfail(reason="To be removed")
def test_doppler_gamma():
    assert False


@pytest.mark.xfail(reason="To be removed")
def test_angle_aberration_gamma():
    assert False


@pytest.mark.parametrize(
    ["energy", "expected"],
    [
        (511.0, 1.0000021334560507),
        (255.5, 0.5000010667280254),
        (0.0, 0.0),
        (511.0e7, 10000021.334560508),
    ],
)
def test_kappa_calculation(energy, expected):
    """

    Parameters
    ----------
    energy : float
    expected : float
    """
    kappa = util.kappa_calculation(energy)
    npt.assert_almost_equal(kappa, expected)


@pytest.mark.xfail(reason="To be removed")
def test_euler_rodrigues():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_solve_quadratic_equation():
    assert False


@pytest.mark.parametrize(
    ["energy", "theta_C"],
    [
        (511.0e3, 1.0),
        (255.5e3, np.pi),
        (0.0, 2.0 * np.pi),
        (511.0e10, np.pi / 2.0),
    ],
)
def test_klein_nishina(energy, theta_C):
    """

    Parameters
    ----------
    energy : float
    theta_C : float
        In radians
    """
    actual = util.klein_nishina(energy, theta_C)

    kappa = util.kappa_calculation(energy)

    expected = (
        R_ELECTRON_SQUARED
        / 2
        * (1.0 + kappa * (1.0 - np.cos(theta_C))) ** -2.0
        * (
            1.0
            + np.cos(theta_C) ** 2.0
            + (kappa ** 2.0 * (1.0 - np.cos(theta_C)) ** 2.0)
            / (1.0 + kappa * (1.0 - np.cos(theta_C)))
        )
    )

    npt.assert_almost_equal(actual, expected)


@pytest.mark.xfail(reason="To be implemented")
def test_compton_theta_distribution():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_normalize_vector():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_get_perpendicular_vector():
    assert False
