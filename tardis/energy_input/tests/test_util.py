import numpy as np
import numpy.testing as npt
import pytest
from numpy.typing import NDArray

from tardis.energy_input.util import (
    C_CGS,
    R_ELECTRON_SQUARED,
    angle_aberration_gamma,
    doppler_factor_3d_all_packets,
    euler_rodrigues,
    get_perpendicular_vector,
    klein_nishina,
    solve_quadratic_equation,
    spherical_to_cartesian,
)
from tardis.opacities.opacities import (
    kappa_calculation,
)


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
    actual_x, actual_y, actual_z = spherical_to_cartesian(r, theta, phi)
    npt.assert_almost_equal(actual_x, expected_x)
    npt.assert_almost_equal(actual_y, expected_y)
    npt.assert_almost_equal(actual_z, expected_z)


def test_doppler_factor_3d_all_packets() -> None:
    """Test Doppler factors for multiple packet vectors."""
    directions = np.array(
        [[1.0, 0.0], [0.0, 0.6], [0.0, 0.8]], dtype=np.float64
    )
    positions = np.array(
        [[3.0e13, 1.0e13], [0.0, 4.0e13], [0.0, 3.0e13]], dtype=np.float64
    )
    times = np.array([1.0e5, 2.0e5], dtype=np.float64)

    expected = 1 - np.sum(positions / times * directions, axis=0) / C_CGS

    npt.assert_allclose(
        doppler_factor_3d_all_packets(directions, positions, times), expected
    )


def test_angle_aberration_gamma_inverse_round_trip() -> None:
    direction = np.array([0.8, 0.6, 0.0])
    position = np.array([3.0e13, 0.0, 0.0])
    time = 2.0e5

    comoving_direction = angle_aberration_gamma(direction, position, time)
    actual = angle_aberration_gamma(
        comoving_direction, position, time, inverse=True
    )
    previous = angle_aberration_gamma(comoving_direction, position, -time)

    npt.assert_allclose(actual, direction)
    npt.assert_allclose(actual, previous)


@pytest.mark.parametrize(
    ("rotation_axis", "initial_direction", "expected_direction"),
    [
        (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ),
        (
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0]),
        ),
        (
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        ),
    ],
)
def test_euler_rodrigues(
    rotation_axis: NDArray[np.float64],
    initial_direction: NDArray[np.float64],
    expected_direction: NDArray[np.float64],
) -> None:
    """Verify right-handed quarter turns about the Cartesian basis axes.

    Claim: A positive pi/2 rotation maps each basis direction according to the
    right-hand rule.
    Regime: Unit rotation axes and unit input vectors in Cartesian space.
    Verification: Exact analytic rotations of the Cartesian basis.
    """
    rotation_matrix = euler_rodrigues(np.pi / 2.0, rotation_axis)

    npt.assert_allclose(
        rotation_matrix @ initial_direction,
        expected_direction,
        rtol=1.0e-14,
        atol=1.0e-15,
        err_msg="Euler-Rodrigues rotation violated the right-hand rule",
    )


@pytest.mark.parametrize(
    ("position", "direction", "radius", "expected_roots"),
    [
        (
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            2.0,
            np.array([1.0, -3.0]),
        ),
        (
            np.array([-3.0, 2.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            2.0,
            np.array([3.0, -np.inf]),
        ),
        (
            np.array([3.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            2.0,
            np.array([-np.inf, -np.inf]),
        ),
    ],
)
def test_solve_quadratic_equation(
    position: NDArray[np.float64],
    direction: NDArray[np.float64],
    radius: float,
    expected_roots: NDArray[np.float64],
) -> None:
    """Verify line-sphere intersections against analytic roots.

    Claim: The solver returns the path lengths where a ray intersects a sphere.
    Regime: A secant, tangent, and non-intersecting unit-direction ray.
    Verification: Closed-form roots of the corresponding quadratic equations.
    """
    actual_roots = solve_quadratic_equation(position, direction, radius)

    npt.assert_allclose(
        actual_roots,
        expected_roots,
        rtol=1.0e-14,
        err_msg="Ray path lengths do not match the analytic sphere intersections",
    )


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
    actual = klein_nishina(energy, theta_C)

    kappa = kappa_calculation(energy)

    expected = (
        R_ELECTRON_SQUARED
        / 2
        * (1.0 + kappa * (1.0 - np.cos(theta_C))) ** -2.0
        * (
            1.0
            + np.cos(theta_C) ** 2.0
            + (kappa**2.0 * (1.0 - np.cos(theta_C)) ** 2.0)
            / (1.0 + kappa * (1.0 - np.cos(theta_C)))
        )
    )

    npt.assert_almost_equal(actual, expected)


def test_get_perpendicular_vector():
    """Test the perpendicular vector calculation"""
    input_vector = np.array([0.3, 0.4, 0.5])

    random_perpendicular_vector = get_perpendicular_vector(input_vector)

    dot_product = np.dot(input_vector, random_perpendicular_vector)

    npt.assert_almost_equal(dot_product, 0.0)
