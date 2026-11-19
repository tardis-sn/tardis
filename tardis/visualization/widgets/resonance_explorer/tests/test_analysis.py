"""Regression and resonance checks for the single-shell analysis."""

import math

import numpy.testing as npt
import pytest

from tardis.conftest import sync_ndarray_assert_allclose
from tardis.visualization.widgets.resonance_explorer import (
    ExplorerInputs,
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.geometry import (
    calculate_projected_gradient,
    calculate_projected_velocity,
)

# Launch point where a path with impact parameter 1.2 enters the sphere r = 3.
OUTER_EDGE_LAUNCH = -math.sqrt(3.0**2 - 1.2**2)
# Turns sit at r³ = -q * b² / m; q = -2.34375 puts them at x = ±0.9 (r = 1.5).
INTERCEPT_FOR_TURNS_AT_0_9 = -2.34375


def make_inputs(
    velocity_slope: float,
    velocity_intercept: float,
    line_velocity: float,
    b: float = 1.2,
    x_0: float = OUTER_EDGE_LAUNCH,
    s_stop: float | None = None,
) -> ExplorerInputs:
    """Build a shell from r=1 to r=3 with v(r) = slope * r + intercept."""
    return ExplorerInputs(
        1.0,
        3.0,
        velocity_slope + velocity_intercept,
        3.0 * velocity_slope + velocity_intercept,
        b,
        x_0,
        line_velocity,
        s_stop,
    )


GENERAL = make_inputs(1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0.3)

EDGE_CASES = [
    pytest.param(GENERAL, id="general"),
    pytest.param(
        make_inputs(1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0.3, b=0, x_0=-3),
        id="radial_path",
    ),
    pytest.param(make_inputs(1.0, 0, 0.3), id="homology"),
    pytest.param(
        make_inputs(1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0), id="line_centre"
    ),
    pytest.param(make_inputs(0, 0.8, 0.3), id="constant_speed"),
    pytest.param(make_inputs(0, 0, 0), id="static_shell"),
    pytest.param(
        make_inputs(1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0.3, b=1, x_0=0),
        id="inner_sphere_tangent_at_launch",
    ),
    pytest.param(
        make_inputs(1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0.7),
        id="target_above_velocity_range",
    ),
    pytest.param(
        make_inputs(
            1.0, INTERCEPT_FOR_TURNS_AT_0_9, 0.275, x_0=-1.6, s_stop=3.2
        ),
        id="scatter_cuts_later_crossing",
    ),
    pytest.param(
        ExplorerInputs(1, 3, 1, 3, 0, 1, 2, s_stop=2),
        id="scatter_exactly_at_boundary",
    ),
]


@pytest.mark.parametrize("inputs", EDGE_CASES)
def test_analysis_matches_regression_data(
    inputs: ExplorerInputs, regression_data: object
) -> None:
    """Claim: distances, crossings and turns are reproducible."""
    case = analyze_case(inputs)
    crossings = case.crossings
    sync_ndarray_assert_allclose(
        regression_data,
        [case.boundary_distance, case.active_end, len(crossings)],
        [crossing.x for crossing in crossings],
        [crossing.s for crossing in crossings],
        [crossing.projected_velocity for crossing in crossings],
        [crossing.gradient for crossing in crossings],
        [turn.x for turn in case.turns],
    )


@pytest.mark.parametrize("inputs", EDGE_CASES)
def test_crossings_satisfy_resonance_condition(inputs: ExplorerInputs) -> None:
    """Claim: every reported crossing has v_parallel = w on the active path."""
    case = analyze_case(inputs)
    for crossing in case.crossings:
        npt.assert_allclose(
            calculate_projected_velocity(crossing.x, inputs),
            inputs.w,
            atol=1e-12,
        )
        assert 0 < crossing.s <= case.active_end
        assert inputs.r_in <= crossing.r <= inputs.r_out


@pytest.mark.parametrize("x", [-2.5, -0.9, 0.0, 0.9, 2.5])
def test_projected_gradient_is_velocity_derivative(x: float) -> None:
    """Claim: G(x) is the derivative of v_parallel(x) along the path."""
    step = 1e-6
    central_difference = (
        calculate_projected_velocity(x + step, GENERAL)
        - calculate_projected_velocity(x - step, GENERAL)
    ) / (2 * step)
    npt.assert_allclose(
        calculate_projected_gradient(x, GENERAL), central_difference, atol=1e-8
    )


@pytest.mark.parametrize(
    "changes",
    [
        dict(r_in=0),
        dict(r_out=0.5),
        dict(b=-1),
        dict(x_0=4),
        dict(s_stop=20),
    ],
)
def test_invalid_geometry_is_actionable(changes: dict[str, float]) -> None:
    """Claim: impossible shells and flights fail before solver evaluation."""
    values = dict(
        r_in=GENERAL.r_in,
        r_out=GENERAL.r_out,
        v_in=GENERAL.v_in,
        v_out=GENERAL.v_out,
        b=GENERAL.b,
        x_0=GENERAL.x_0,
        w=GENERAL.w,
    )
    values.update(changes)
    with pytest.raises(ValueError, match=r"(?i)shell|launch|stop|impact"):
        analyze_case(ExplorerInputs(**values))
