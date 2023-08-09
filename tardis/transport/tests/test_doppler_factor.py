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
    # Set the params from test cases here
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = frame_transformations.get_doppler_factor(r, mu, time_explosion)

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
    # Set the params from test cases here
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = frame_transformations.get_inverse_doppler_factor(
        r, mu, time_explosion
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
    obtained = frame_transformations.get_inverse_doppler_factor_full_relativity(
        mu, beta
    )
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(["mu", "r", "inv_t_exp"], [(-0.3, 5, 1e10)])
def test_unphysical_doppler_factor(mu, r, inv_t_exp):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables
    with pytest.raises(r_packet.SuperluminalError):
        obtained = frame_transformations.get_doppler_factor(
            r, mu, time_explosion
        )


@pytest.mark.parametrize(["mu", "r", "inv_t_exp"], [(-0.3, 5, 1e10)])
def test_unphysical_inverse_doppler_factor(mu, r, inv_t_exp):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables
    with pytest.raises(r_packet.SuperluminalError):
        obtained = r_packet.get_inverse_doppler_factor(r, mu, time_explosion)
