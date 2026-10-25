"""Reconstruct structured-solver branch predicates beside its actual slots.

The production solver does not expose a runtime trace. These values reproduce
its source expressions and invoke its private depressed-quartic routine only
for direct-root comparison. They are explicitly a reconstructed diagnostic.
"""

import math

from tardis.transport.montecarlo.nonhomologous_grid import (
    FLOAT64_ROUNDOFF_TOLERANCE,
    FLOAT64_TINY,
    RECIPROCAL_ROOT_THRESHOLD,
    _solve_depressed_quartic,
    depressed_quartic,
)

from .analysis import CaseAnalysis


def diagnose_solver(analysis: CaseAnalysis) -> dict[str, object]:
    """Return source-matched branch decisions for one analysis snapshot.

    Parameters
    ----------
    analysis : CaseAnalysis
        Immutable analysis containing actual production root slots.

    Returns
    -------
    dict[str, object]
        Labeled reconstructed predicates and comparison with actual slots.
    """
    base: dict[str, object] = {
        "provenance": "reconstructed diagnostic; actual slots retained separately",
        "route": analysis.solver_route,
        "actual_slots": analysis.raw_roots,
        "factorization": "not entered",
        "discriminant_sign": "not evaluated",
        "reciprocal_triggered": False,
        "reciprocal_attempted": False,
        "reciprocal_selected": False,
        "agrees_with_solver": None,
        "comparison_status": "not reconstructed for this top-level route",
    }
    if analysis.solver_route != "general quartic":
        return base

    structured = dict(analysis.structured_parameters)
    line_velocity = structured["N"]
    impact_squared = structured["H²"]
    intercept = structured["Q"]
    variable_scale = max(
        abs(line_velocity), math.sqrt(impact_squared), abs(intercept)
    )
    scaled_line = line_velocity / variable_scale
    scaled_impact_squared = impact_squared / (variable_scale * variable_scale)
    scaled_intercept = intercept / variable_scale

    line_squared = scaled_line * scaled_line
    intercept_squared = scaled_intercept * scaled_intercept
    depressed_quadratic = (
        scaled_impact_squared - intercept_squared - 0.5 * line_squared
    )
    depressed_linear = -scaled_line * (
        scaled_impact_squared + intercept_squared
    )
    depressed_constant = (
        line_squared
        * (4.0 * (scaled_impact_squared - intercept_squared) + line_squared)
        / 16.0
    )
    invariant_sum = scaled_impact_squared + line_squared - intercept_squared
    invariant_product = scaled_impact_squared * line_squared * intercept_squared
    inner = invariant_sum**3 + 27.0 * invariant_product
    inner_scale = abs(invariant_sum**3) + 27.0 * invariant_product
    clamped = abs(inner) <= FLOAT64_ROUNDOFF_TOLERANCE * max(
        inner_scale, FLOAT64_TINY
    )
    if clamped:
        inner = 0.0
    discriminant = 4.0 * invariant_product * inner / 27.0
    discriminant_sign = (
        "positive"
        if discriminant > 0
        else "negative"
        if discriminant < 0
        else "zero"
    )

    direct_roots = _solve_depressed_quartic(
        depressed_quadratic,
        depressed_linear,
        depressed_constant,
        0.5 * scaled_line,
        discriminant,
        True,
    )
    real_count = sum(math.isfinite(root) for root in direct_roots)
    reciprocal_triggered = any(
        math.isfinite(root) and abs(root) < RECIPROCAL_ROOT_THRESHOLD
        for root in direct_roots
    )
    monic_quadratic = line_squared + scaled_impact_squared - intercept_squared
    monic_linear = -2.0 * scaled_line * scaled_impact_squared
    monic_constant = line_squared * scaled_impact_squared
    reciprocal_attempted = reciprocal_triggered and monic_constant != 0.0
    reciprocal_roots: tuple[float, float, float, float] | None = None
    selected_roots = tuple(float(root) for root in direct_roots)
    reciprocal_selected = False
    if reciprocal_attempted:
        reciprocal_roots = depressed_quartic(
            monic_constant,
            monic_linear,
            monic_quadratic,
            -2.0 * scaled_line,
            1.0,
        )
        mixed = [
            root
            for root in direct_roots
            if math.isfinite(root) and abs(root) >= RECIPROCAL_ROOT_THRESHOLD
        ]
        for reciprocal_root in reciprocal_roots:
            if math.isfinite(reciprocal_root) and reciprocal_root != 0.0:
                root = 1.0 / reciprocal_root
                if abs(root) < RECIPROCAL_ROOT_THRESHOLD and len(mixed) < 4:
                    mixed.append(root)
        if len(mixed) == real_count:
            selected_roots = tuple((mixed + [math.nan] * (4 - len(mixed)))[:4])
            reciprocal_selected = True

    reconstructed = tuple(variable_scale * root for root in selected_roots)
    slots_agree = all(
        (math.isnan(first) and math.isnan(second))
        or math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-14)
        for first, second in zip(reconstructed, analysis.raw_roots, strict=True)
    )
    base.update(
        {
            "variable_scale": variable_scale,
            "depressed_quadratic": depressed_quadratic,
            "depressed_linear": depressed_linear,
            "depressed_constant": depressed_constant,
            "factorization": "biquadratic"
            if depressed_linear == 0.0
            else "Ferrari",
            "discriminant_inner_clamped": clamped,
            "resolvent_discriminant": discriminant,
            "discriminant_sign": discriminant_sign,
            "direct_scaled_slots": tuple(float(root) for root in direct_roots),
            "reciprocal_triggered": reciprocal_triggered,
            "reciprocal_attempted": reciprocal_attempted,
            "reciprocal_scaled_slots": reciprocal_roots,
            "reciprocal_selected": reciprocal_selected,
            "reconstructed_slots": reconstructed,
            "agrees_with_solver": slots_agree,
            "comparison_status": "ordered reconstructed slots compared with actual solver slots",
        }
    )
    return base


GENERIC_SOLVER_EXAMPLES: dict[str, tuple[float, float, float, float, float]] = {
    "biquadratic four real": (1.0, 0.0, -5.0, 0.0, 4.0),
    "positive resolvent discriminant": (1.0, 0.0, 0.0, 1.0, 0.0),
    "generic zero real": (1.0, 0.0, 0.0, 1.0, 1.0),
    "repeated real contact": (1.0, 0.0, -1.0, -2.0, 2.0),
}


def diagnose_generic_quartic(
    coefficients: tuple[float, float, float, float, float],
) -> dict[str, object]:
    """Summarize a solver-only quartic outside the physical input model.

    Parameters
    ----------
    coefficients : tuple[float, float, float, float, float]
        Quartic coefficients from highest to lowest power, with nonzero lead.

    Returns
    -------
    dict[str, object]
        Generic factorization route, discriminant branch, and actual slots.
    """
    leading, cubic, quadratic, linear, constant = coefficients
    if not all(math.isfinite(value) for value in coefficients) or leading == 0:
        raise ValueError(
            "Generic diagnostic requires finite quartic coefficients and nonzero lead"
        )
    monic_cubic = cubic / leading
    monic_quadratic = quadratic / leading
    monic_linear = linear / leading
    monic_constant = constant / leading
    p = monic_quadratic - 3.0 * monic_cubic**2 / 8.0
    q = (
        monic_cubic**3 / 8.0
        - monic_cubic * monic_quadratic / 2.0
        + monic_linear
    )
    r = (
        -3.0 * monic_cubic**4 / 256.0
        + monic_quadratic * monic_cubic**2 / 16.0
        - monic_cubic * monic_linear / 4.0
        + monic_constant
    )
    discriminant: float | None = None
    sign = "not evaluated"
    if q != 0.0:
        alpha = -(p * p + 12.0 * r) / 3.0
        beta = -(2.0 * p**3 - 72.0 * p * r + 27.0 * q * q) / 27.0
        positive_term = beta * beta / 4.0
        negative_term = alpha**3 / 27.0
        discriminant = positive_term + negative_term
        tolerance = FLOAT64_ROUNDOFF_TOLERANCE * max(
            abs(positive_term), abs(negative_term), FLOAT64_TINY
        )
        if abs(discriminant) <= tolerance:
            discriminant = 0.0
        sign = (
            "positive"
            if discriminant > 0
            else "negative"
            if discriminant < 0
            else "zero"
        )
    slots = depressed_quartic(*coefficients)
    return {
        "provenance": "solver-only generic quartic",
        "coefficients": coefficients,
        "depressed_coefficients": (p, q, r),
        "factorization": "biquadratic" if q == 0 else "Ferrari",
        "resolvent_discriminant": discriminant,
        "discriminant_sign": sign,
        "actual_slots": slots,
        "finite_root_count": sum(math.isfinite(root) for root in slots),
    }
