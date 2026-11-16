"""Verification for source-matched, reconstructed branch predicates."""

from math import sqrt

from tardis.visualization.widgets.resonance_explorer.analysis import (
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.diagnostics import (
    GENERIC_SOLVER_EXAMPLES,
    diagnose_generic_quartic,
    diagnose_solver,
)
from tardis.visualization.widgets.resonance_explorer.model import ExplorerInputs
from tardis.visualization.widgets.resonance_explorer.presets import PRESETS


def test_linked_case_uses_negative_discriminant_ferrari_path() -> None:
    """Claim: the linked four-slot case enters the negative-resolvent branch."""
    snapshot = analyze_case(PRESETS["linked three crossings"].inputs)
    diagnostic = diagnose_solver(snapshot)
    assert diagnostic["factorization"] == "Ferrari"
    assert diagnostic["discriminant_sign"] == "negative"
    assert diagnostic["agrees_with_solver"] is True


def test_zero_intercept_skips_general_factorization() -> None:
    """Claim: the exact Q=0 route precedes general quartic scaling."""
    snapshot = analyze_case(PRESETS["radius proportional"].inputs)
    diagnostic = diagnose_solver(snapshot)
    assert diagnostic["factorization"] == "not entered"
    assert diagnostic["discriminant_sign"] == "not evaluated"


def test_generic_zero_real_example_is_solver_only() -> None:
    """Claim: a zero-real generic quartic is not a physical shell preset."""
    diagnostic = diagnose_generic_quartic(
        GENERIC_SOLVER_EXAMPLES["generic zero real"]
    )
    assert diagnostic["provenance"] == "solver-only generic quartic"
    assert diagnostic["finite_root_count"] == 0


def test_generic_biquadratic_has_four_real_slots() -> None:
    """Claim: the generic biquadratic branch returns ±1 and ±2."""
    diagnostic = diagnose_generic_quartic(
        GENERIC_SOLVER_EXAMPLES["biquadratic four real"]
    )
    assert diagnostic["factorization"] == "biquadratic"
    assert diagnostic["finite_root_count"] == 4


def test_generic_examples_expose_each_discriminant_sign() -> None:
    """Claim: solver-only examples span positive, zero, and negative signs."""
    signs = {
        diagnose_generic_quartic(coefficients)["discriminant_sign"]
        for coefficients in GENERIC_SOLVER_EXAMPLES.values()
    }
    assert {"positive", "zero", "negative"} <= signs


def test_reciprocal_recovery_is_reported_for_tiny_physical_roots() -> None:
    """Claim: source threshold triggers reciprocal recovery for tiny roots."""
    impact = 0.00233329371
    intercept = -973.326786
    inputs = ExplorerInputs(
        1.0,
        2.0,
        1.0 + intercept,
        2.0 + intercept,
        impact,
        -sqrt(4.0 - impact**2),
        -0.0113652295,
    )
    diagnostic = diagnose_solver(analyze_case(inputs))
    assert diagnostic["reciprocal_triggered"] is True
    assert diagnostic["reciprocal_attempted"] is True
    assert diagnostic["reciprocal_selected"] is True
    assert diagnostic["agrees_with_solver"] is True
