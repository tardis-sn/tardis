import numpy as np
import numpy.testing as ntest
import pytest


from tardis.spectrum.formal_integral.base import check, intensity_black_body
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.spectrum.formal_integral.formal_integral_numba import (
    calculate_p_values as calculate_p_values_numba,
    intensity_black_body as intensity_black_body_numba,
)
from tardis.spectrum.formal_integral.formal_integral_cuda import (
    calculate_p_values as calculate_p_values_cuda,
    intensity_black_body_cuda,
)


@pytest.mark.parametrize(
    "line_interaction_type",
    ("downbranch", "macroatom", pytest.param("?", marks=pytest.mark.xfail)),
)
def test_check(simulation_verysimple, line_interaction_type):
    sim_state = simulation_verysimple.simulation_state
    plasma = simulation_verysimple.plasma
    transport = simulation_verysimple.transport
    transport.line_interaction_type = line_interaction_type

    assert check(sim_state, plasma, transport)

    # should return false
    assert not check(None, plasma, transport, raises=False)
    assert not check(sim_state, None, transport, raises=False)
    assert not check(sim_state, plasma, None, raises=False)


@pytest.mark.parametrize(
    ["nu", "temperature", "expected"],
    [
        (10**6, 1000, 3.072357852080765e-22),
        (10**6, 300, 9.21707305730458e-23),
        (10**8, 1000, 6.1562660718558254e-24),
        (10**8, 300, 1.846869480674048e-24),
    ],
)
def test_intensity_black_body(nu, temperature, expected):
    actual = intensity_black_body(nu, temperature)
    assert np.isclose(actual, expected)

    actual_numba = intensity_black_body_numba(nu, temperature)
    assert np.isclose(actual_numba, expected)

    # TODO: check if cuda
    # actual_cuda = intensity_black_body_cuda(nu, temperature)
    # assert np.isclose(actual_cuda, expected)


@pytest.mark.parametrize(
    "N",
    [
        100,
        1000,
        10000,
    ],
)
def test_calculate_p_values(N):
    r = 1.0

    expected = r / (N - 1) * np.arange(0, N, dtype=np.float64)
    actual = np.zeros_like(expected, dtype=np.float64)

    actual[::] = calculate_p_values_numba(r, N)
    ntest.assert_allclose(actual, expected)

    # TODO: check if cuda
    # actual_cuda = np.zeros_like(expected, dtype=np.float64)
    # actual_cuda[::] = calculate_p_values_cuda(r, N)
    # ntest.assert_allclose(actual_cuda, expected)
