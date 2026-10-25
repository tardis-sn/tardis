"""Physics and navigation contracts for the seven guided lessons."""

import math
from itertools import pairwise

import numpy.testing as npt
from scipy.optimize import brentq

from tardis.visualization.widgets.resonance_explorer.analysis import (
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.diagnostics import (
    diagnose_solver,
)
from tardis.visualization.widgets.resonance_explorer.lessons import (
    LESSONS,
    get_lesson,
    get_lesson_state,
)
from tardis.visualization.widgets.resonance_explorer.presets import PRESETS


def test_guided_navigation_has_seven_distinct_questions() -> None:
    """Claim: a newcomer can enter the seven physics questions directly."""
    assert tuple(lesson.key for lesson in LESSONS) == (
        "homology",
        "linked",
        "launch",
        "target",
        "collapse",
        "homologous_collapse",
        "branches",
    )
    assert all(lesson.question and lesson.states for lesson in LESSONS)
    assert all(get_lesson(lesson.key) is lesson for lesson in LESSONS)


def test_homology_comparison_changes_only_radial_intercept() -> None:
    """Claim: lessons 1 and 2 isolate the intercept at fixed shell and ray."""
    homology = get_lesson_state("homology", "initial").inputs
    linked = get_lesson_state("linked", "initial").inputs
    npt.assert_allclose(
        (homology.r_in, homology.r_out, homology.b, homology.x_0, homology.w),
        (linked.r_in, linked.r_out, linked.b, linked.x_0, linked.w),
        rtol=0,
        atol=0,
    )
    npt.assert_allclose(homology.x_0, -math.sqrt(7.56), atol=1e-15)
    npt.assert_allclose((homology.m, linked.m), (1.0, 1.0), atol=0)
    npt.assert_allclose((homology.q, linked.q), (0.0, -2.34375), atol=0)
    npt.assert_allclose(
        (homology.v_in, homology.v_out, linked.v_in, linked.v_out),
        (1.0, 3.0, -1.34375, 0.65625),
        atol=0,
    )


def test_homologous_lessons_have_one_analytic_crossing_and_constant_gradient() -> (
    None
):
    """Claim: q=0 gives v_parallel=m*x and G=m on the full flyby."""
    for lesson_key, m, w in (
        ("homology", 1.0, 0.3),
        ("homologous_collapse", -1.0, 0.6),
    ):
        inputs = get_lesson_state(lesson_key, "initial").inputs
        analysis = analyze_case(inputs)
        assert analysis.solver_route == "zero-intercept Q=0"
        assert len(analysis.turns) == 0
        assert len(analysis.crossings) == 1
        npt.assert_allclose(analysis.crossings[0].x, w / m, atol=1e-14)
        npt.assert_allclose(analysis.crossings[0].gradient, m, atol=1e-14)
        finite_slots = [
            root for root in analysis.raw_roots if math.isfinite(root)
        ]
        assert len(finite_slots) == 2  # one repeated algebraic root
        npt.assert_allclose(finite_slots, [w / (m * inputs.delta_r)] * 2)


def test_launch_and_target_states_change_one_causal_input() -> None:
    """Claim: lessons 3 and 4 hold the linked velocity law and ray fixed."""
    linked = get_lesson_state("linked", "initial").inputs
    for state in get_lesson("launch").states:
        assert state.highlight_input == "x_0"
        assert state.inputs.x_0 != linked.x_0
        assert state.inputs.b == linked.b
        assert state.inputs.w == linked.w
        assert state.inputs.v_in == linked.v_in
        assert state.inputs.v_out == linked.v_out
    for state in get_lesson("target").states:
        assert state.highlight_input == "w"
        assert state.inputs.w != linked.w
        assert state.inputs.b == linked.b
        assert state.inputs.x_0 == linked.x_0
        assert state.inputs.v_in == linked.v_in
        assert state.inputs.v_out == linked.v_out


def test_linked_launch_and_target_crossing_counts() -> None:
    """Claim: clipping and moving w produce the intended active topology."""
    expected = (
        ("linked", "initial", 2, 3),
        ("launch", "two_crossings", 1, 2),
        ("launch", "one_crossing", 0, 1),
        ("target", "above_range", 2, 0),
        ("target", "touch", 2, 2),
    )
    for lesson_key, state_key, turn_count, crossing_count in expected:
        analysis = analyze_case(get_lesson_state(lesson_key, state_key).inputs)
        assert analysis.solver_route == "general quartic"
        assert len(analysis.turns) == turn_count
        assert len(analysis.crossings) == crossing_count
    linked = analyze_case(get_lesson_state("linked", "initial").inputs)
    assert len(linked.raw_roots) == 4
    assert (
        sum("squared-sign artifact" in c.flags for c in linked.candidates) == 1
    )
    touch = analyze_case(get_lesson_state("target", "touch").inputs)
    npt.assert_allclose(touch.analytic_contacts[0].x, -0.9, atol=1e-12)
    touch_inputs = get_lesson_state("target", "touch").inputs
    contact_radius = math.hypot(-0.9, touch_inputs.b)
    npt.assert_allclose(
        -0.9 + touch_inputs.q * -0.9 / contact_radius,
        touch_inputs.w,
        atol=1e-14,
    )
    npt.assert_allclose(
        1 + touch_inputs.q * touch_inputs.b**2 / contact_radius**3,
        0,
        atol=1e-14,
    )


def test_wholly_collapsing_shell_has_two_independently_bracketed_crossings() -> (
    None
):
    """Claim: negative radial speed can coexist with two inbound resonances.

    Regime: m=1, q=-3.5, b=1.2, w=0.6 on the full shell flyby.
    Verification: bracket the original unsquared equation on its monotone
    intervals, independently of the production quartic and squared equation.
    """
    inputs = get_lesson_state("collapse", "initial").inputs
    analysis = analyze_case(inputs)
    assert inputs.m == 1.0
    assert inputs.q == -3.5
    assert inputs.v_in == -2.5
    assert inputs.v_out == -0.5
    assert analysis.solver_route == "general quartic"

    def original_equation(x: float) -> float:
        return x + (-3.5) * x / math.hypot(x, 1.2) - 0.6

    turn = math.sqrt((3.5 * 1.2**2) ** (2.0 / 3.0) - 1.2**2)
    boundaries = (inputs.x_0, -turn, turn, -inputs.x_0)
    bracketed = [
        brentq(original_equation, left, right, xtol=1e-13)
        for left, right in pairwise(boundaries)
        if original_equation(left) * original_equation(right) < 0
    ]
    # Eight displayed decimal places limit agreement with the brief's
    # rounded roots to about 5e-9; Brent's numerical root is tighter.
    npt.assert_allclose(
        bracketed, (-2.57170109, -0.33120996), atol=6e-9, rtol=0
    )
    npt.assert_allclose(
        [candidate.x for candidate in analysis.crossings],
        bracketed,
        atol=1e-9,
        rtol=0,
    )
    assert len(analysis.turns) == 2
    assert all(candidate.x < 0 for candidate in analysis.crossings)

    outside = next(
        candidate
        for candidate in analysis.candidates
        if candidate.x is not None and candidate.x > inputs.x_0 * -1
    )
    assert "outside shell" in outside.flags
    outside_root = brentq(original_equation, 3.8, 4.1, xtol=1e-13)
    npt.assert_allclose(outside_root, 3.94878442, atol=6e-9, rtol=0)
    npt.assert_allclose(outside.x, outside_root, atol=1e-9, rtol=0)
    artifact = next(
        candidate
        for candidate in analysis.candidates
        if "squared-sign artifact" in candidate.flags
    )
    npt.assert_allclose(artifact.x, 0.15412663, atol=6e-9, rtol=0)
    assert abs(original_equation(artifact.x)) > 0.8
    squared_equation = (artifact.x - inputs.w) ** 2 * (
        artifact.x**2 + inputs.b**2
    ) - inputs.q**2 * artifact.x**2
    npt.assert_allclose(squared_equation, 0, atol=1e-12)


def test_branch_sampler_reuses_reference_cases_and_positive_discriminant() -> (
    None
):
    """Claim: exact bypass routes and a physical positive branch are present."""
    for state_key, preset_name, expected_route in (
        ("radial", "radial inward", "radial H²=0"),
        ("uniform", "uniform speed analytic", "Δv=0 analytic"),
        ("line_center", "line center", "line-center N=0"),
        ("scatter", "scatter truncation", "general quartic"),
    ):
        state = get_lesson_state("branches", state_key)
        assert state.inputs == PRESETS[preset_name].inputs
        assert analyze_case(state.inputs).solver_route == expected_route
    positive = get_lesson_state("branches", "positive_discriminant")
    assert positive.inputs.m == 1.0
    assert positive.inputs.q == 0.5
    assert positive.inputs.w == 0.3
    positive_analysis = analyze_case(positive.inputs)
    assert positive_analysis.solver_route == "general quartic"
    assert len(positive_analysis.crossings) == 1
    assert positive.inputs.m > 0 and positive.inputs.q > 0

    # G = m + q*b²/r³ is positive everywhere on this ray, so this
    # independently bracketed zero is the only possible physical crossing.
    def positive_equation(x: float) -> float:
        return x + 0.5 * x / math.hypot(x, 1.2) - 0.3

    unique_root = brentq(
        positive_equation, positive.inputs.x_0, -positive.inputs.x_0
    )
    npt.assert_allclose(
        positive_analysis.crossings[0].x, unique_root, atol=1e-12
    )
    assert diagnose_solver(positive_analysis)["discriminant_sign"] == "positive"
    linked = analyze_case(get_lesson_state("linked", "initial").inputs)
    assert diagnose_solver(linked)["discriminant_sign"] == "negative"
