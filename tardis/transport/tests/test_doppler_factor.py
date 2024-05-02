import pytest
import numpy as np

import tardis.transport.frame_transformations as frame_transformations
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet

from numpy.testing import (
    assert_almost_equal,
)


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (0.3, 7.5e14, 1 / 5.2e7, 0.9998556693818854),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_doppler_factor(mu, r, inv_t_exp, expected):
    """
    Checks the get_doppler_factor function.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    r : float
        Radius of the position of the packet.
    inv_t_exp : float
        Inverse of t_explosion.
    expected : float
        Expected value of the doppler factor.
    """
    # Set the params from test cases here
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = frame_transformations.get_doppler_factor(
        r, mu, time_explosion, False
    )

    # Perform required assertions
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "beta", "expected"],
    [
        (0.3, 0.2, 0.94),
        (-0.3, 0, 1.0),
        (0, 0.8, 1.0),
    ],
)
def test_get_doppler_factor_partial_relativity(mu, beta, expected):
    """
    Checks the get_doppler_factor_partial_relativity.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    beta : float
        Velocity over speed of light for the packet.
    expected : float
        Expected value of the doppler factor.
    """
    obtained = frame_transformations.get_doppler_factor_partial_relativity(
        mu, beta
    )
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "beta", "expected"],
    [
        (0.3, 0.2, 0.95938348),
        (-0.3, 0, 1.0),
        (0, 0.8, 1.6666667),
    ],
)
def test_get_doppler_factor_full_relativity(mu, beta, expected):
    """
    Checks the get_doppler_factor_full_relativity.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    beta : float
        Velocity over speed of light for the packet.
    expected : float
        Expected value of the doppler factor.
    """
    obtained = frame_transformations.get_doppler_factor_full_relativity(
        mu, beta
    )
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (0.3, 7.5e14, 1 / 5.2e7, 1 / 0.9998556693818854),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_inverse_doppler_factor(mu, r, inv_t_exp, expected):
    """
    Checks the get_inverse_doppler_factor function.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    r : float
        Radius of the position of the packet.
    inv_t_exp : float
        Inverse of t_explosion.
    expected : float
        Expected value of the inverse doppler factor.
    """
    # Set the params from test cases here
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = frame_transformations.get_inverse_doppler_factor(
        r, mu, time_explosion, enable_full_relativity=False
    )

    # Perform required assertions
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "beta", "expected"],
    [
        (0.3, 0.2, 1 / 0.94),
        (-0.3, 0, 1.0),
        (0, 0.8, 1.0),
    ],
)
def test_get_inverse_doppler_factor_partial_relativity(mu, beta, expected):
    """
    Checks the get_inverse_doppler_factor_partial_relativity function.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    beta : float
        Velocity over speed of light for the packet.
    expected : float
        Expected value of the inverse doppler factor.
    """
    obtained = (
        frame_transformations.get_inverse_doppler_factor_partial_relativity(
            mu, beta
        )
    )
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "beta", "expected"],
    [
        (0.3, 0.2, 1.0818579),
        (-0.3, 0, 1.0),
        (0, 0.8, 1.6666667),
    ],
)
def test_get_inverse_doppler_factor_full_relativity(mu, beta, expected):
    """
    Checks the get_inverse_doppler_factor_full_relativity function.

    Parameters
    ----------
    mu : float
        Angle of movement of the packet.
    beta : float
        Velocity over speed of light for the packet.
    expected : float
        Expected value of the inverse doppler factor.
    """
    obtained = frame_transformations.get_inverse_doppler_factor_full_relativity(
        mu, beta
    )
    assert_almost_equal(obtained, expected)
