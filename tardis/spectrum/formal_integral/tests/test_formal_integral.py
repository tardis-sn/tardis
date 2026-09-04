from types import SimpleNamespace

import numpy as np
import numpy.testing as ntest
import pandas as pd
import pytest

from tardis.spectrum.formal_integral.base import (
    check_formal_integral_requirements,
    intensity_black_body,
)
from tardis.spectrum.formal_integral.formal_integral_numba import (
    calculate_impact_parameters as calculate_impact_parameters_numba,
)
from tardis.spectrum.formal_integral.formal_integral_numba import (
    intensity_black_body as intensity_black_body_numba,
)
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)
from tardis.spectrum.formal_integral.source_function import SourceFunctionState


@pytest.mark.parametrize(
    "line_interaction_type",
    ("downbranch", "macroatom", pytest.param("?", marks=pytest.mark.xfail)),
)
def test_check_formal_integral_requirements(
    simulation_verysimple, line_interaction_type
):
    sim_state = simulation_verysimple.simulation_state
    plasma = simulation_verysimple.plasma
    transport = simulation_verysimple.transport
    transport.line_interaction_type = line_interaction_type

    assert check_formal_integral_requirements(sim_state, plasma, transport)

    # should return false
    warning_match = (
        "The integrator is missing either model, opacity state or transport"
    )
    with pytest.warns(UserWarning, match=warning_match):
        assert not check_formal_integral_requirements(
            None, plasma, transport, raises=False
        )
    with pytest.warns(UserWarning, match=warning_match):
        assert not check_formal_integral_requirements(
            sim_state, None, transport, raises=False
        )
    with pytest.warns(UserWarning, match=warning_match):
        assert not check_formal_integral_requirements(
            sim_state, plasma, None, raises=False
        )


def test_interpolate_integrator_quantities_uses_active_shells() -> None:
    solver = FormalIntegralSolver(points=10, interpolate_shells=-1)
    source_function_state = SourceFunctionState(
        att_S_ul=np.ones((2, 2)),
        Jred_lu=np.ones((2, 2)),
        Jblue_lu=np.ones((2, 2)),
        e_dot_u=pd.DataFrame(),
    )
    opacity_state = SimpleNamespace(tau_sobolev=pd.DataFrame(np.ones((2, 2))))

    result = solver.interpolate_integrator_quantities(
        r_inner_original=np.array([1.0, 2.0]),
        r_outer_original=np.array([2.0, 3.0]),
        r_inner_interpolated=np.array([1.0, 1.5, 2.0]),
        r_outer_interpolated=np.array([1.5, 2.0, 2.5]),
        source_function_state=source_function_state,
        opacity_state=opacity_state,
        electron_densities=pd.Series([10.0, 20.0]),
    )

    assert result[5].shape == (2, 3)
    assert result[6].shape == (3,)


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

    actual[::] = calculate_impact_parameters_numba(r, N)
    ntest.assert_allclose(actual, expected)

    # TODO: check if cuda
    # actual_cuda = np.zeros_like(expected, dtype=np.float64)
    # actual_cuda[::] = calculate_impact_parameters_cuda(r, N)
    # ntest.assert_allclose(actual_cuda, expected)
