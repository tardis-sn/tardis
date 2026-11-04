"""Independent equation checks for representative named manuscript cases.

These tests use the original unsquared equation and fixed geometric brackets.
They do not infer expected topology from preset metadata, quartic slots, or
the explorer's computed monotonic intervals.
"""

import math
from dataclasses import replace

import numpy as np

from tardis.visualization.widgets.resonance_explorer import analyze_case
from tardis.visualization.widgets.resonance_explorer.presets import PRESETS


def _bisect_unsquared_root(lower: float, upper: float, target: float) -> float:
    """Bisect the manuscript's m=1, q=-2.34375, b=1.2 equation."""

    def residual(x: float) -> float:
        return x - 2.34375 * x / math.hypot(x, 1.2) - target

    lower_value = residual(lower)
    assert lower_value * residual(upper) < 0
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if lower_value * residual(midpoint) <= 0:
            upper = midpoint
        else:
            lower = midpoint
            lower_value = residual(lower)
    return (lower + upper) / 2


def test_linked_preset_matches_three_independent_unsquared_brackets() -> None:
    """Claim: the published full-chord example has three ordered roots.

    Regime: m=1, q=-2.34375, b=1.2, w=0.3; no earlier stop.
    Verification: bisect the original equation in three fixed disjoint
    brackets, then compare to manuscript coordinates rounded to six places.
    The 1e-6 manuscript tolerance reflects that printed precision; the
    solver-to-bisection tolerance is 1e-12 length units.
    """
    case = analyze_case(PRESETS["linked three crossings"].inputs)
    bracketed = np.asarray(
        [
            _bisect_unsquared_root(-1.8, -1.3, 0.3),
            _bisect_unsquared_root(-0.5, -0.1, 0.3),
            _bisect_unsquared_root(2.2, 2.5, 0.3),
        ]
    )
    np.testing.assert_allclose(
        bracketed, [-1.555887, -0.341473, 2.395531], atol=1e-6, rtol=0
    )
    assert len(case.crossings) == 3
    np.testing.assert_allclose(
        [crossing.x for crossing in case.crossings],
        bracketed,
        atol=1e-12,
        rtol=0,
    )


def test_named_target_miss_exceeds_direct_active_extrema() -> None:
    """Claim: w=0.7 is unreachable on the linked full chord.

    Verification: G=1-2.34375*b²/(x²+b²)^(3/2) vanishes at x=±0.9.
    Compare the target to those two extrema and both shell endpoints.
    """
    case = analyze_case(PRESETS["target above range"].inputs)
    edge_x = math.sqrt(3**2 - 1.2**2)
    values = [
        x - 2.34375 * x / math.hypot(x, 1.2)
        for x in (-edge_x, -0.9, 0.9, edge_x)
    ]
    assert case.inputs.w > max(values)
    np.testing.assert_allclose(max(values), 0.601463059962954, atol=1e-12)
    assert case.crossings == ()
    np.testing.assert_allclose(
        case.active_velocity_range, (min(values), max(values))
    )


def test_named_extremum_touch_has_one_contact_and_one_simple_root() -> None:
    """Claim: the crest target touches at -0.9 and cuts once outward.

    Verification: direct substitution gives f(-0.9)=0 and G(-0.9)=0;
    independent bisection locates the simple right root. The 1e-7 contact
    tolerance covers the production solver's repeated-root splitting.
    """
    case = analyze_case(PRESETS["extremum touch"].inputs)
    contact_x = -0.9
    target = 0.50625
    np.testing.assert_allclose(
        contact_x - 2.34375 * contact_x / math.hypot(contact_x, 1.2),
        target,
        atol=0,
        rtol=0,
    )
    np.testing.assert_allclose(
        1 - 2.34375 * 1.2**2 / math.hypot(contact_x, 1.2) ** 3,
        0,
        atol=1e-15,
    )
    right_root = _bisect_unsquared_root(2.6, 2.7, target)
    np.testing.assert_allclose(right_root, 2.639909, atol=1e-6)
    assert len(case.crossings) == 2
    np.testing.assert_allclose(case.crossings[0].x, contact_x, atol=1e-7)
    np.testing.assert_allclose(case.crossings[1].x, right_root, atol=1e-12)


def test_close_simple_roots_remain_distinct_below_crest() -> None:
    """Claim: a target 1e-5 below the crest cuts both sides of the peak.

    Verification: independent brackets straddle x=-0.9 and enclose two
    roots separated by about 0.008 length units, far above roundoff.
    """
    inputs = PRESETS["linked three crossings"].inputs
    case = analyze_case(replace(inputs, w=0.50624))
    left_root = _bisect_unsquared_root(-0.91, -0.9, 0.50624)
    right_root = _bisect_unsquared_root(-0.9, -0.89, 0.50624)
    assert 0.005 < right_root - left_root < 0.01
    assert len(case.crossings) == 3
    np.testing.assert_allclose(
        [case.crossings[0].x, case.crossings[1].x],
        [left_root, right_root],
        atol=1e-10,
        rtol=0,
    )


def test_named_scatter_excludes_later_independently_bracketed_root() -> None:
    """Claim: an assumed scatter at x=1.6 clips the right resonance.

    Verification: direct bisection gives one root before and one after the
    stated stop; the latter remains a verified raw slot but cannot be a
    numbered active crossing.
    """
    case = analyze_case(PRESETS["scatter truncation"].inputs)
    before = _bisect_unsquared_root(-0.5, -0.1, 0.275)
    after = _bisect_unsquared_root(2.2, 2.5, 0.275)
    np.testing.assert_allclose(case.x_end, 1.6, atol=1e-15)
    assert before < case.x_end < after
    assert len(case.crossings) == 1
    np.testing.assert_allclose(case.crossings[0].x, before, atol=1e-12)
    assert any(
        candidate.physical
        and candidate.x is not None
        and abs(candidate.x - after) < 1e-12
        and candidate.s > case.active_end
        for candidate in case.candidates
    )


def test_named_radial_inward_root_follows_unsquared_linear_limit() -> None:
    """Claim: an inward radial packet resonates at x=-2.04375.

    Verification: b=0 and x<0 reduce the unsquared equation to
    v_parallel=x+2.34375; w=0.3 gives x=w-2.34375 exactly.
    """
    case = analyze_case(PRESETS["radial inward"].inputs)
    expected_x = 0.3 - 2.34375
    assert case.inputs.x_0 < expected_x < -case.inputs.r_in
    assert len(case.crossings) == 1
    np.testing.assert_allclose(case.crossings[0].x, expected_x, atol=1e-14)
    np.testing.assert_allclose(case.crossings[0].residual, 0, atol=1e-14)
