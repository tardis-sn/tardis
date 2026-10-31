"""Packet-facing arithmetic checks for structured solver inputs."""

import math

import numpy as np
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

from ..analysis import _classify_candidate  # noqa: TID252
from ..model import ExplorerInputs  # noqa: TID252 - Package-local test import.
from ..solver_adapter import solve_candidates  # noqa: TID252


@pytest.mark.parametrize(
    ("r_in", "r_out", "v_in", "v_out", "expected_route"),
    [
        (0.2, 1.3, 0.13999999999999999, 0.9099999999999999, "general quartic"),
        (
            1.1,
            3.2,
            0.11000000000000001,
            0.32000000000000006,
            "zero-intercept Q=0",
        ),
    ],
)
def test_exact_q_zero_route_matches_packet_arithmetic(
    r_in: float,
    r_out: float,
    v_in: float,
    v_out: float,
    expected_route: str,
) -> None:
    """Claim: the displayed intercept and Q=0 route use packet arithmetic.

    Regime: decimal edge values where inner-edge and outer-edge subtractions
    differ by one floating-point unit. Verification: use the production
    geometry's gradient and the packet distance routine's exact outer-edge
    intercept expression; exact zero is a branch predicate, not a tolerance.
    """
    geometry = NumbaRadial1DGeometry(
        np.array([r_in]),
        np.array([r_out]),
        np.array([v_in]),
        np.array([v_out]),
    )
    packet_q = (
        geometry.v_outer[0]
        - geometry.velocity_gradient[0] * geometry.r_outer[0]
    )
    impact = 0.1
    inputs = ExplorerInputs(
        r_in,
        r_out,
        v_in,
        v_out,
        impact,
        -math.sqrt(r_out**2 - impact**2),
        0.2,
    )
    solver = solve_candidates(inputs)

    assert inputs.q == packet_q
    assert dict(solver.structured_parameters)["Q"] == packet_q / (v_out - v_in)
    assert solver.route == expected_route


def test_near_radial_h_squared_uses_packet_direction_arithmetic() -> None:
    """Claim: a near-radial packet takes the exact production H²=0 route.

    Regime: b=1e-9 at x_0=2 makes the stored packet cosine exactly one.
    Verification: independently reconstruct the scaled impact parameter from
    radius and cosine, as calculate_distance_line_nonhomologous does.
    The exact zero comparison matters because it selects a solver branch.
    """
    inputs = ExplorerInputs(1, 3, -1.34375, 0.65625, 1e-9, 2, 0.3)
    packet_radius = math.hypot(inputs.x_0, inputs.b)
    packet_cosine = inputs.x_0 / packet_radius
    packet_h_squared = (packet_radius / inputs.delta_r) ** 2 * (
        1.0 - packet_cosine**2
    )
    solver = solve_candidates(inputs)

    assert packet_h_squared == 0.0
    assert dict(solver.structured_parameters)["H²"] == packet_h_squared
    assert solver.route == "radial H²=0"


def test_shell_filter_uses_packet_position_arithmetic_at_narrow_edge() -> None:
    """Claim: shell membership follows reconstructed packet radius.

    Regime: a narrow shell whose geometric outer-edge `hypot` rounds to the
    boundary while the packet-position expression rounds one ulp beyond it.
    Verification: evaluate the source's radius expression independently;
    the shell tolerance is much smaller than one ulp in this case.
    """
    r_outer = 1.00001
    b = 0.8
    x_0 = math.sqrt(1.0 - b * b)
    x_edge = math.sqrt(r_outer * r_outer - b * b)
    inputs = ExplorerInputs(1.0, r_outer, 1.0, r_outer, b, x_0, x_edge)
    s = x_edge - x_0
    packet_radius = math.sqrt(
        inputs.r_0 * inputs.r_0 + s * s + 2 * inputs.r_0 * s * inputs.mu_0
    )
    tolerance = 1e-14 * inputs.delta_r

    assert math.hypot(x_edge, b) == r_outer
    assert packet_radius > r_outer + tolerance
    candidate = _classify_candidate(
        0, x_edge / inputs.delta_r, inputs, 1.0, 1.0, tolerance
    )
    assert "outside shell" in candidate.flags
    assert not candidate.accepted
