"""Packet-facing arithmetic checks for structured solver inputs."""

import math

import numpy as np
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.visualization.widgets.resonance_explorer.model import ExplorerInputs
from tardis.visualization.widgets.resonance_explorer.solver_adapter import (
    solve_candidates,
)


@pytest.mark.parametrize(
    ("r_inner", "r_outer", "v_inner", "v_outer", "expected_route"),
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
    r_inner: float,
    r_outer: float,
    v_inner: float,
    v_outer: float,
    expected_route: str,
) -> None:
    """Claim: the displayed intercept and Q=0 route use packet arithmetic.

    Regime: decimal edge values where inner-edge and outer-edge subtractions
    differ by one floating-point unit. Verification: use the production
    geometry's gradient and the packet distance routine's exact outer-edge
    intercept expression; exact zero is a branch predicate, not a tolerance.
    """
    geometry = NumbaRadial1DGeometry(
        np.array([r_inner]),
        np.array([r_outer]),
        np.array([v_inner]),
        np.array([v_outer]),
    )
    packet_q = (
        geometry.v_outer[0]
        - geometry.velocity_gradient[0] * geometry.r_outer[0]
    )
    impact = 0.1
    inputs = ExplorerInputs(
        r_inner,
        r_outer,
        v_inner,
        v_outer,
        impact,
        -math.sqrt(r_outer**2 - impact**2),
        0.2,
    )
    solver = solve_candidates(inputs)

    assert inputs.velocity_intercept == packet_q
    assert solver.route == expected_route


def test_near_radial_h_squared_uses_packet_direction_arithmetic() -> None:
    """Claim: a near-radial packet takes the exact production H²=0 route.

    Regime: impact_parameter=1e-9 at launch_position=2 makes the stored packet cosine exactly one.
    Verification: independently reconstruct the scaled impact parameter from
    radius and cosine, as calculate_distance_line_nonhomologous does.
    The exact zero comparison matters because it selects a solver branch.
    """
    inputs = ExplorerInputs(1, 3, -1.34375, 0.65625, 1e-9, 2, 0.3)
    packet_radius = math.hypot(inputs.launch_position, inputs.impact_parameter)
    packet_cosine = inputs.launch_position / packet_radius
    packet_h_squared = (packet_radius / inputs.delta_r) ** 2 * (
        1.0 - packet_cosine**2
    )
    solver = solve_candidates(inputs)

    assert packet_h_squared == 0.0
    assert solver.route == "radial H²=0"
