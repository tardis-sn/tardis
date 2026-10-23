"""Scientific checks for the single-shell resonance analysis."""

import math

import numpy as np
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_line_nonhomologous,
    calculate_projected_gradient_zero_distances,
)
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    CLOSE_LINE_THRESHOLD,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket
from tardis.visualization.widgets.resonance_explorer import (
    ExplorerInputs,
    analyze_case,
)


def linked_inputs(**changes: float) -> ExplorerInputs:
    """Return the paper's linked three-crossing geometry."""
    values = dict(
        r_in=1.0,
        r_out=3.0,
        v_in=-1.34375,
        v_out=0.65625,
        b=1.2,
        x_0=-math.sqrt(3.0**2 - 1.2**2),
        w=0.3,
    )
    values.update(changes)
    return ExplorerInputs(**values)


def test_linked_case_has_three_unsquared_crossings_and_two_turns() -> None:
    """Claim: one full flyby can cross one target three times.

    Verification: the published roots and analytic G=0 positions are
    independent of the production quartic's slot order.
    """
    case = analyze_case(linked_inputs())
    np.testing.assert_allclose(
        sorted(root.x for root in case.crossings),
        [-1.555887, -0.341473, 2.395531],
        atol=1e-6,
    )
    np.testing.assert_allclose([turn.x for turn in case.turns], [-0.9, 0.9])
    assert case.ray_family == "inner-miss flyby"
    assert case.shape_family == "two turns"
    assert case.solver_route == "general quartic"
    assert len(case.raw_roots) == 4


def test_linked_crossings_match_independent_interval_bisection() -> None:
    """Claim: each selected root solves the original unsquared equation.

    Verification: bracket one root in each analytically monotone segment
    and bisect the direct projected-velocity equation without the quartic.
    """
    case = analyze_case(linked_inputs())

    def direct_residual(position: float) -> float:
        return position - 2.34375 * position / math.hypot(position, 1.2) - 0.3

    bracketed = []
    for lower, upper in case.monotonic_intervals:
        lower_value = direct_residual(lower)
        upper_value = direct_residual(upper)
        assert lower_value * upper_value < 0
        for _ in range(60):
            midpoint = (lower + upper) / 2
            if direct_residual(midpoint) * lower_value <= 0:
                upper = midpoint
            else:
                lower = midpoint
                lower_value = direct_residual(lower)
        bracketed.append((lower + upper) / 2)

    np.testing.assert_allclose(
        [root.x for root in case.crossings], bracketed, atol=1e-12
    )


def test_ill_conditioned_slot_reports_unsquared_verification() -> None:
    """Claim: production acceptance alone does not certify a line crossing.

    Regime: a nearly uniform, fast inward radial law makes a structured
    quartic slot inaccurate even though its production sign filter passes.
    Verification: bisect the original equation on its monotone flight and
    compare positions using the local derivative. The forward error budget
    is 64 machine epsilons times the absolute terms in that equation; its
    position equivalent is that budget divided by the local gradient.
    """
    inputs = ExplorerInputs(
        1,
        3,
        -1103.1408330241957,
        -1103.137596530609,
        1.0987592278970888,
        -2.7915458368278663,
        22.32587042642611,
    )
    case = analyze_case(inputs)
    candidate = case.candidates[1]

    def unsquared_residual(x: float) -> float:
        return inputs.m * x + inputs.q * x / math.hypot(x, inputs.b) - inputs.w

    # Fixed bracket from the direct equation, independent of explorer splits.
    lower, upper = -0.03, -0.01
    assert inputs.x_0 < lower < upper < case.x_end
    assert unsquared_residual(lower) * unsquared_residual(upper) < 0
    lower_residual = unsquared_residual(lower)
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if unsquared_residual(midpoint) * lower_residual <= 0:
            upper = midpoint
        else:
            lower = midpoint
            lower_residual = unsquared_residual(lower)
    bracketed_x = (lower + upper) / 2
    term_scale = (
        abs(inputs.m * candidate.x)
        + abs(inputs.q * candidate.x / math.hypot(candidate.x, inputs.b))
        + abs(inputs.w)
    )
    forward_error = 64 * np.finfo(float).eps * term_scale
    position_tolerance = forward_error / abs(candidate.gradient)
    assert abs(candidate.x - bracketed_x) > 1000 * position_tolerance
    assert candidate.accepted
    assert candidate.event_eligible
    assert case.production_nearest_resonance == candidate
    assert case.crossings == ()
    assert case.nearest_resonance is None
    assert not candidate.unsquared_verified
    assert not candidate.verified_event_eligible


def test_nearest_resonance_matches_packet_distance_at_physical_scale() -> None:
    """Claim: nearest-root selection matches the packet-facing distance API.

    Verification: map lengths to 10¹³ cm and velocities to 10⁸ cm/s so
    the line-frequency ratio remains representable, then call TARDIS.
    """
    case = analyze_case(linked_inputs())
    length_unit = 1.0e13
    velocity_unit = 1.0e8
    packet_frequency = 1.0e15
    geometry = NumbaRadial1DGeometry(
        np.asarray([case.inputs.r_in * length_unit]),
        np.asarray([case.inputs.r_out * length_unit]),
        np.asarray([case.inputs.v_in * velocity_unit]),
        np.asarray([case.inputs.v_out * velocity_unit]),
    )
    packet = RPacket(
        r=case.r_0 * length_unit,
        mu=case.mu_0,
        nu=packet_frequency,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    line_frequency = packet_frequency * (
        1 - case.inputs.w * velocity_unit / C_SPEED_OF_LIGHT
    )
    packet_distance = calculate_distance_line_nonhomologous(
        packet,
        geometry,
        line_frequency,
        maximum_distance=case.active_end * length_unit,
    )
    np.testing.assert_allclose(
        packet_distance, case.nearest_resonance.s * length_unit, rtol=1e-10
    )


def test_radial_inward_root_has_no_gradient_reversal() -> None:
    """Claim: radial inward flow has constant projected gradient m.

    Verification: x=w+q is the exact inward-leg root for m=1.
    """
    case = analyze_case(ExplorerInputs(1, 3, -1.34375, 0.65625, 0, -3, 0.3))
    assert case.boundary_kind == "inner"
    assert case.ray_family == "inner-boundary approach"
    assert case.solver_route == "radial H²=0"
    assert len(case.turns) == 0
    np.testing.assert_allclose([root.x for root in case.crossings], [-2.04375])


def test_constant_speed_analytic_and_static_nonisolated() -> None:
    """Claim: a uniform-speed nonradial ray has one analytic isolated root.

    Verification: squaring q*x/sqrt(x²+b²)=w gives the root magnitude,
    and the original sign determines its sign.
    """
    case = analyze_case(ExplorerInputs(1, 3, 2, 2, 1.2, -2, 1))
    np.testing.assert_allclose(
        [root.x for root in case.crossings], [1.2 / math.sqrt(3)]
    )
    assert case.solver_route == "Δv=0 analytic"
    static_case = analyze_case(ExplorerInputs(1, 3, 0, 0, 1.2, -2, 0))
    assert static_case.result_family == "non-isolated match"
    assert len(static_case.crossings) == 0


def test_negative_uniform_radial_speed_match_is_nonisolated() -> None:
    """Claim: the matched constant-speed radial target follows flow sign."""
    outward = analyze_case(ExplorerInputs(1, 3, -2, -2, 0, 1, -2))
    inward = analyze_case(ExplorerInputs(1, 3, -2, -2, 0, -3, 2))
    assert outward.result_family == "non-isolated match"
    assert inward.result_family == "non-isolated match"


def test_exact_inner_tangent_uses_production_boundary_ownership() -> None:
    """Claim: production boundary code owns an exact inner tangent at launch."""
    case = analyze_case(ExplorerInputs(1, 3, 1, 3, 1, 0, 0))
    assert case.boundary_kind == "inner"
    assert case.boundary_distance == 0
    assert case.x_boundary == 0
    assert case.active_end == 0
    assert case.ray_family == "outbound leg"
    assert case.shape_family == "flat"
    assert case.monotonic_intervals == ()
    assert any("geometric outbound continuation" in note for note in case.notes)


def test_scatter_truncates_later_physical_root() -> None:
    """Claim: an assumed earlier scatter truncates the active flight."""
    case = analyze_case(linked_inputs(s_stop=2.0))
    assert case.active_end == 2.0
    assert all(root.s <= 2 for root in case.crossings)
    assert any(
        "after event stop" in candidate.flags for candidate in case.candidates
    )


def test_target_touch_is_one_distinct_singular_contact() -> None:
    """Claim: the crest target has one repeated root and one simple root.

    Verification: G vanishes analytically at x=-0.9, while the other
    unsquared root is the paper's x≈2.639909. Raw slots retain multiplicity.
    """
    case = analyze_case(linked_inputs(w=0.50625))
    assert len(case.crossings) == 2
    assert len(case.analytic_contacts) == 1
    np.testing.assert_allclose(case.analytic_contacts[0].x, -0.9)
    np.testing.assert_allclose(case.crossings[-1].x, 2.639909, atol=1e-6)
    assert case.result_family == "extremum touch plus simple crossing"


def test_adjacent_target_is_not_labeled_exact_contact() -> None:
    """Claim: a representable near-contact target has distinct provenance."""
    near_target = math.nextafter(0.50625, math.inf)
    case = analyze_case(linked_inputs(w=near_target))
    assert case.analytic_contacts == ()
    assert not case.result_family.startswith("extremum touch")


def test_flat_center_crossing_is_not_a_turn() -> None:
    """Claim: q=-m*b gives a stationary inflection at x=0.

    Verification: the analytic gradient has one zero with no sign change;
    the fourfold polynomial slot represents one unsquared crossing.
    """
    case = analyze_case(linked_inputs(v_in=-0.2, v_out=1.8, w=0))
    assert case.shape_family == "stationary inflection"
    assert case.result_family == "flat-center crossing"
    assert len(case.turns) == 0
    assert len(case.crossings) == 1
    # The analytic flat-center limit is exact, but packet-form arithmetic
    # retains a finite G for these decimal edge speeds.
    assert case.crossings[0].packet_gradient != 0.0
    assert "singular projected gradient" not in case.crossings[0].flags
    assert any("packet-form G=" in note for note in case.notes)


def test_small_nonzero_gradient_is_not_singular() -> None:
    """Claim: packet transport's singular-depth condition is exactly G==0.

    Regime: proportional radial flow with an eligible root and G≈1e-13.
    Verification: G=m=1e-13 analytically for b=0 and q=0; the packet
    transport branch compares the computed gradient to zero exactly.
    """
    case = analyze_case(ExplorerInputs(1, 3, 1e-13, 3e-13, 0, 1, 2e-13))
    assert len(case.crossings) == 1
    assert case.crossings[0].gradient != 0.0
    assert "singular projected gradient" not in case.crossings[0].flags


def test_exact_zero_packet_gradient_is_singular() -> None:
    """Claim: an exact G=0 isolated crossing has singular Sobolev depth.

    Verification: m=1, q=-b=-2 gives a stationary inflection at x=0;
    integer edge speeds make packet v(r=2)=0 exactly, so its gradient is
    exactly zero in the same arithmetic used by transport.
    """
    case = analyze_case(ExplorerInputs(1, 3, -1, 1, 2, -math.sqrt(5), 0))
    assert len(case.crossings) == 1
    assert case.crossings[0].packet_gradient == 0.0
    assert "singular projected gradient" in case.crossings[0].flags


@pytest.mark.parametrize(
    ("epsilon_multiple", "expected_transport_turns"),
    [(0.25, 0), (2.5, 2)],
)
def test_near_endpoint_turn_has_distinct_geometric_and_transport_status(
    epsilon_multiple: float, expected_transport_turns: int
) -> None:
    """Claim: a geometric turn splits transport only beyond epsilon_d.

    Verification: compare the analytic G=0 construction with the packet
    routine on both sides of epsilon_d=CLOSE_LINE_THRESHOLD*shell width.
    The separation is larger than floating-point rounding at radius 3.
    """
    b = 1.2
    x_0 = -math.sqrt(3**2 - b**2)
    epsilon_d = CLOSE_LINE_THRESHOLD * (3 - 1)
    negative_turn_x = x_0 + epsilon_multiple * epsilon_d
    q = -((negative_turn_x**2 + b**2) ** 1.5) / b**2
    inputs = ExplorerInputs(1, 3, 1 + q, 3 + q, b, x_0, 0)
    case = analyze_case(inputs)
    geometry = NumbaRadial1DGeometry(
        np.asarray([inputs.r_in], dtype=float),
        np.asarray([inputs.r_out], dtype=float),
        np.asarray([inputs.v_in], dtype=float),
        np.asarray([inputs.v_out], dtype=float),
    )
    packet = RPacket(
        r=case.r_0,
        mu=case.mu_0,
        nu=1e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    _, _, production_count = calculate_projected_gradient_zero_distances(
        packet, geometry, case.boundary_distance
    )
    assert len(case.turns) == 2
    assert production_count == expected_transport_turns
    assert len(case.transport_turns) == production_count
    assert len(case.monotonic_intervals) == production_count + 1
    assert all(
        turn.transport_active == (expected_transport_turns == 2)
        for turn in case.turns
    )


def test_boundary_tie_and_scatter_tie_have_distinct_owners() -> None:
    """Claim: the shell owns an exact boundary tie; a scatter tie leaves the line.

    Verification: proportional radial flow has v_parallel=x on the
    outward leg, making both event distances exact integers.
    """
    at_boundary = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, 1, 3))
    assert len(at_boundary.crossings) == 1
    assert at_boundary.nearest_resonance is None
    assert any(
        "boundary tie" in candidate.flags
        for candidate in at_boundary.candidates
    )
    at_scatter = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, 1, 2, s_stop=1))
    assert at_scatter.nearest_resonance is not None
    np.testing.assert_allclose(at_scatter.nearest_resonance.s, 1)


def test_target_at_launch_is_excluded_by_lower_distance_tolerance() -> None:
    """Claim: a launch coincidence is not a forward line event."""
    case = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, 1, 1))
    assert len(case.crossings) == 0
    assert any(
        "at or behind start" in candidate.flags for candidate in case.candidates
    )


def test_no_crossing_reports_target_outside_active_velocity_range() -> None:
    """Claim: a target above every active value explains a missed line.

    Verification: compare w=0.7 to the outer-boundary value 0.601463,
    which exceeds the interior crest 0.50625.
    """
    case = analyze_case(linked_inputs(w=0.7))
    assert case.crossings == ()
    np.testing.assert_allclose(case.active_velocity_range[1], 0.601463059962954)
    assert "outside active projected-velocity range" in case.no_crossing_reason


def test_no_crossing_explains_verified_roots_outside_active_flight() -> None:
    """Claim: an earlier scatter can leave a real root beyond the active ray.

    Verification: for proportional radial flow v_parallel=x, w=2.5 occurs
    at x=2.5, while the packet is stopped at x=2.
    """
    case = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, 1, 2.5, s_stop=1))
    assert case.crossings == ()
    assert any(
        candidate.physical and candidate.s > case.active_end
        for candidate in case.candidates
    )
    assert "verified root" in case.no_crossing_reason.lower()
    assert "outside active flight" in case.no_crossing_reason


def test_scatter_at_boundary_does_not_claim_truncation() -> None:
    """Claim: an exact scatter/boundary tie retains boundary ownership.

    Verification: radial outward flight from r=1 to r=3 has boundary
    distance 2, identical to the assumed scatter stop.
    """
    case = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, 1, 2, s_stop=2))
    assert case.stop_reason == "outer"
    assert not any("truncated by" in note for note in case.notes)


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
    with pytest.raises(ValueError, match=r"(?i)shell|launch|stop|impact"):
        analyze_case(linked_inputs(**changes))
