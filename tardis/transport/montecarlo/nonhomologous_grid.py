import math

import numpy as np
from numba import njit

_FLOAT64_EPSILON = np.finfo(np.float64).eps
_FLOAT64_TINY = np.finfo(np.float64).tiny
_RECIPROCAL_ROOT_THRESHOLD = math.sqrt(_FLOAT64_EPSILON)
_QUARTIC_NJIT_OPTIONS = {
    "fastmath": False,
    "error_model": "numpy",
    "parallel": False,
}


@njit(**_QUARTIC_NJIT_OPTIONS)
def _solve_real_quadratic(
    linear_coefficient: float, constant_coefficient: float
) -> tuple[float, float]:
    """Solve a monic quadratic without cancellation between its roots."""
    discriminant = (
        linear_coefficient * linear_coefficient - 4.0 * constant_coefficient
    )
    discriminant_scale = max(
        abs(linear_coefficient * linear_coefficient),
        abs(4.0 * constant_coefficient),
        _FLOAT64_TINY,
    )
    discriminant_tolerance = 8.0 * _FLOAT64_EPSILON * discriminant_scale

    if discriminant < 0.0:
        if discriminant >= -discriminant_tolerance:
            discriminant = 0.0
        else:
            return math.nan, math.nan

    square_root = math.sqrt(discriminant)
    if square_root == 0.0:
        repeated_root = -0.5 * linear_coefficient
        return repeated_root, repeated_root

    first_root = -0.5 * (
        linear_coefficient + math.copysign(square_root, linear_coefficient)
    )
    if first_root == 0.0:
        return (
            0.5 * (-linear_coefficient + square_root),
            0.5 * (-linear_coefficient - square_root),
        )
    return first_root, constant_coefficient / first_root


@njit(**_QUARTIC_NJIT_OPTIONS)
def _solve_biquadratic(
    quadratic_coefficient: float,
    constant_coefficient: float,
    root_shift: float,
) -> tuple[float, float, float, float]:
    """Solve a depressed quartic with no linear term."""
    first_squared_root, second_squared_root = _solve_real_quadratic(
        quadratic_coefficient, constant_coefficient
    )
    squared_root_tolerance = (
        8.0
        * _FLOAT64_EPSILON
        * max(
            abs(quadratic_coefficient),
            math.sqrt(abs(constant_coefficient)),
            _FLOAT64_TINY,
        )
    )
    first_positive_root = math.nan
    first_negative_root = math.nan
    second_positive_root = math.nan
    second_negative_root = math.nan

    if first_squared_root < 0.0:
        if first_squared_root >= -squared_root_tolerance:
            first_squared_root = 0.0
        else:
            first_positive_root = math.nan
            first_negative_root = math.nan
    if first_squared_root >= 0.0:
        square_root = math.sqrt(first_squared_root)
        first_positive_root = root_shift + square_root
        first_negative_root = root_shift - square_root

    if second_squared_root < 0.0:
        if second_squared_root >= -squared_root_tolerance:
            second_squared_root = 0.0
        else:
            second_positive_root = math.nan
            second_negative_root = math.nan
    if second_squared_root >= 0.0:
        square_root = math.sqrt(second_squared_root)
        second_positive_root = root_shift + square_root
        second_negative_root = root_shift - square_root

    return (
        first_positive_root,
        first_negative_root,
        second_positive_root,
        second_negative_root,
    )


@njit(**_QUARTIC_NJIT_OPTIONS)
def _solve_resolvent_cubic(
    quadratic_coefficient: float,
    linear_coefficient: float,
    constant_coefficient: float,
    discriminant: float,
    use_supplied_discriminant: bool,
) -> float:
    """Return the largest real solution for ``W**2`` in the factorization."""
    cubic_quadratic_coefficient = 2.0 * quadratic_coefficient
    depressed_linear_coefficient = (
        -(
            quadratic_coefficient * quadratic_coefficient
            + 12.0 * constant_coefficient
        )
        / 3.0
    )
    depressed_constant_coefficient = (
        -(
            2.0 * quadratic_coefficient**3
            - 72.0 * quadratic_coefficient * constant_coefficient
            + 27.0 * linear_coefficient * linear_coefficient
        )
        / 27.0
    )

    if not use_supplied_discriminant:
        positive_discriminant_term = (
            depressed_constant_coefficient
            * depressed_constant_coefficient
            / 4.0
        )
        negative_discriminant_term = depressed_linear_coefficient**3 / 27.0
        discriminant = positive_discriminant_term + negative_discriminant_term
        discriminant_tolerance = (
            8.0
            * _FLOAT64_EPSILON
            * max(
                abs(positive_discriminant_term),
                abs(negative_discriminant_term),
                _FLOAT64_TINY,
            )
        )
        if abs(discriminant) <= discriminant_tolerance:
            discriminant = 0.0

    cubic_shift = cubic_quadratic_coefficient / 3.0
    if discriminant > 0.0:
        square_root = math.sqrt(discriminant)
        half_constant = -0.5 * depressed_constant_coefficient
        cube_argument = half_constant + math.copysign(
            square_root,
            half_constant if half_constant != 0.0 else 1.0,
        )
        first_cube_root = np.cbrt(cube_argument)
        if first_cube_root == 0.0:
            depressed_root = np.cbrt(-depressed_constant_coefficient)
        else:
            second_cube_root = -depressed_linear_coefficient / (
                3.0 * first_cube_root
            )
            depressed_root = first_cube_root + second_cube_root
        squared_factor = depressed_root - cubic_shift
    elif discriminant == 0.0:
        cube_root = np.cbrt(-0.5 * depressed_constant_coefficient)
        squared_factor = max(
            2.0 * cube_root - cubic_shift,
            -cube_root - cubic_shift,
        )
    else:
        if depressed_linear_coefficient >= 0.0:
            return math.nan
        root_amplitude = 2.0 * math.sqrt(-depressed_linear_coefficient / 3.0)
        cosine_argument = (
            3.0
            * depressed_constant_coefficient
            / (2.0 * depressed_linear_coefficient)
            * math.sqrt(-3.0 / depressed_linear_coefficient)
        )
        cosine_argument = min(1.0, max(-1.0, cosine_argument))
        angle = math.acos(cosine_argument) / 3.0
        squared_factor = -math.inf
        for root_index in range(3):
            candidate = (
                root_amplitude
                * math.cos(angle - 2.0 * math.pi * root_index / 3.0)
                - cubic_shift
            )
            squared_factor = max(squared_factor, candidate)

    squared_factor_scale = max(
        abs(quadratic_coefficient),
        math.sqrt(abs(constant_coefficient)),
        abs(linear_coefficient) ** (2.0 / 3.0),
        _FLOAT64_TINY,
    )
    squared_factor_tolerance = 8.0 * _FLOAT64_EPSILON * squared_factor_scale
    if squared_factor < 0.0:
        if squared_factor >= -squared_factor_tolerance:
            squared_factor = 0.0
        else:
            return math.nan

    # When W**2 is small, the subtraction of the cubic shift above loses
    # relative precision. Recover the same root from the unshifted resolvent
    # identity without subtracting nearly equal terms.
    if squared_factor < 0.01 * squared_factor_scale:
        denominator = (
            quadratic_coefficient + squared_factor
        ) ** 2 - 4.0 * constant_coefficient
        if denominator > 0.0:
            recovered_squared_factor = (
                linear_coefficient * linear_coefficient / denominator
            )
            if recovered_squared_factor > 0.0:
                squared_factor = recovered_squared_factor

    return squared_factor


@njit(**_QUARTIC_NJIT_OPTIONS)
def _solve_depressed_quartic(
    quadratic_coefficient: float,
    linear_coefficient: float,
    constant_coefficient: float,
    root_shift: float,
    resolvent_discriminant: float,
    use_supplied_discriminant: bool,
) -> tuple[float, float, float, float]:
    """Factor and solve a depressed quartic using real arithmetic."""
    if linear_coefficient == 0.0:
        return _solve_biquadratic(
            quadratic_coefficient, constant_coefficient, root_shift
        )

    squared_factor = _solve_resolvent_cubic(
        quadratic_coefficient,
        linear_coefficient,
        constant_coefficient,
        resolvent_discriminant,
        use_supplied_discriminant,
    )
    if not math.isfinite(squared_factor) or squared_factor <= 0.0:
        return math.nan, math.nan, math.nan, math.nan

    factor = math.sqrt(squared_factor)
    factor_sum = quadratic_coefficient + squared_factor
    factor_sum_scale = max(
        abs(quadratic_coefficient), squared_factor, _FLOAT64_TINY
    )
    if abs(factor_sum) < 0.01 * factor_sum_scale:
        factor_sum_squared = (
            linear_coefficient * linear_coefficient / squared_factor
            + 4.0 * constant_coefficient
        )
        factor_sum_squared_tolerance = (
            8.0
            * _FLOAT64_EPSILON
            * max(
                linear_coefficient * linear_coefficient / squared_factor,
                abs(4.0 * constant_coefficient),
                _FLOAT64_TINY,
            )
        )
        if (
            factor_sum_squared < 0.0
            and factor_sum_squared >= -factor_sum_squared_tolerance
        ):
            factor_sum_squared = 0.0
        if factor_sum_squared >= 0.0:
            factor_sum = math.copysign(
                math.sqrt(factor_sum_squared), factor_sum
            )
            squared_factor = factor_sum - quadratic_coefficient
            factor = math.sqrt(squared_factor)
    factor_difference = linear_coefficient / factor
    first_constant = 0.5 * (factor_sum - factor_difference)
    second_constant = 0.5 * (factor_sum + factor_difference)

    # One of these expressions may subtract nearly equal terms. Preserve the
    # larger factor and recover the smaller one from their exact product.
    if abs(first_constant) >= abs(second_constant):
        if first_constant != 0.0:
            second_constant = constant_coefficient / first_constant
    elif second_constant != 0.0:
        first_constant = constant_coefficient / second_constant

    first_root, second_root = _solve_real_quadratic(factor, first_constant)
    third_root, fourth_root = _solve_real_quadratic(-factor, second_constant)
    return (
        root_shift + first_root,
        root_shift + second_root,
        root_shift + third_root,
        root_shift + fourth_root,
    )


@njit(**_QUARTIC_NJIT_OPTIONS)
def depressed_quartic(
    leading_coefficient: float,
    cubic_coefficient: float,
    quadratic_coefficient: float,
    linear_coefficient: float,
    constant_coefficient: float,
) -> tuple[float, float, float, float]:
    """Return the real roots of a quartic and NaN for complex roots.

    The polynomial is normalized before it is depressed, then factored into
    two quadratics using a real solution of its resolvent cubic.

    Parameters
    ----------
    leading_coefficient : float
        Coefficient of the fourth-degree term.
    cubic_coefficient : float
        Coefficient of the third-degree term.
    quadratic_coefficient : float
        Coefficient of the second-degree term.
    linear_coefficient : float
        Coefficient of the first-degree term.
    constant_coefficient : float
        Constant term.

    Returns
    -------
    tuple of float
        Four real roots or NaN placeholders for non-real roots.
    """
    monic_cubic = cubic_coefficient / leading_coefficient
    monic_quadratic = quadratic_coefficient / leading_coefficient
    monic_linear = linear_coefficient / leading_coefficient
    monic_constant = constant_coefficient / leading_coefficient

    depressed_quadratic = monic_quadratic - 3.0 * monic_cubic**2 / 8.0
    depressed_linear = (
        monic_cubic**3 / 8.0
        - monic_cubic * monic_quadratic / 2.0
        + monic_linear
    )
    depressed_constant = (
        -3.0 * monic_cubic**4 / 256.0
        + monic_quadratic * monic_cubic**2 / 16.0
        - monic_cubic * monic_linear / 4.0
        + monic_constant
    )
    root_shift = -monic_cubic / 4.0

    return _solve_depressed_quartic(
        depressed_quadratic,
        depressed_linear,
        depressed_constant,
        root_shift,
        0.0,
        False,
    )


@njit(**_QUARTIC_NJIT_OPTIONS)
def solve_resonance_quartic(
    line_velocity: float,
    impact_parameter_squared: float,
    velocity_intercept: float,
) -> tuple[float, float, float, float]:
    """Solve the structured quartic for a line-resonance distance.

    The polynomial is

    ``(x - N)**2 * (x**2 + H**2) - Q**2 * x**2``,

    where the arguments are ``N``, ``H**2``, and ``Q`` respectively.

    Parameters
    ----------
    line_velocity : float
        Dimensionless line velocity ``N``.
    impact_parameter_squared : float
        Squared dimensionless ray impact parameter ``H**2``.
    velocity_intercept : float
        Dimensionless affine velocity intercept ``Q``.

    Returns
    -------
    tuple of float
        Four real roots or NaN placeholders for non-real roots.
    """
    if impact_parameter_squared == 0.0:
        return (
            line_velocity + velocity_intercept,
            line_velocity - velocity_intercept,
            0.0,
            0.0,
        )

    if velocity_intercept == 0.0:
        return line_velocity, line_velocity, math.nan, math.nan

    if line_velocity == 0.0:
        additional_squared_root = (
            velocity_intercept * velocity_intercept - impact_parameter_squared
        )
        if additional_squared_root < 0.0:
            return 0.0, 0.0, math.nan, math.nan
        additional_root = math.sqrt(additional_squared_root)
        return additional_root, -additional_root, 0.0, 0.0

    variable_scale = max(
        abs(line_velocity),
        math.sqrt(impact_parameter_squared),
        abs(velocity_intercept),
    )
    scaled_line_velocity = line_velocity / variable_scale
    scaled_impact_parameter_squared = impact_parameter_squared / (
        variable_scale * variable_scale
    )
    scaled_velocity_intercept = velocity_intercept / variable_scale

    line_velocity_squared = scaled_line_velocity * scaled_line_velocity
    velocity_intercept_squared = (
        scaled_velocity_intercept * scaled_velocity_intercept
    )
    depressed_quadratic = (
        scaled_impact_parameter_squared
        - velocity_intercept_squared
        - 0.5 * line_velocity_squared
    )
    depressed_linear = -scaled_line_velocity * (
        scaled_impact_parameter_squared + velocity_intercept_squared
    )
    depressed_constant = (
        line_velocity_squared
        * (
            4.0 * (scaled_impact_parameter_squared - velocity_intercept_squared)
            + line_velocity_squared
        )
        / 16.0
    )

    invariant_sum = (
        scaled_impact_parameter_squared
        + line_velocity_squared
        - velocity_intercept_squared
    )
    invariant_product = (
        scaled_impact_parameter_squared
        * line_velocity_squared
        * velocity_intercept_squared
    )
    discriminant_inner_term = invariant_sum**3 + 27.0 * invariant_product
    discriminant_inner_scale = abs(invariant_sum**3) + 27.0 * invariant_product
    if abs(discriminant_inner_term) <= (
        8.0 * _FLOAT64_EPSILON * max(discriminant_inner_scale, _FLOAT64_TINY)
    ):
        discriminant_inner_term = 0.0
    resolvent_discriminant = (
        4.0 * invariant_product * discriminant_inner_term / 27.0
    )

    first_root, second_root, third_root, fourth_root = _solve_depressed_quartic(
        depressed_quadratic,
        depressed_linear,
        depressed_constant,
        0.5 * scaled_line_velocity,
        resolvent_discriminant,
        True,
    )

    direct_roots = (first_root, second_root, third_root, fourth_root)
    recover_small_roots = False
    direct_real_root_count = 0
    for root in direct_roots:
        if math.isfinite(root):
            direct_real_root_count += 1
            if abs(root) < _RECIPROCAL_ROOT_THRESHOLD:
                recover_small_roots = True

    # Shifting depressed roots back to x can cancel when a real root is much
    # smaller than every characteristic scale. In that regime the same root
    # is large and well conditioned in the reciprocal polynomial. Combine
    # large direct roots with small reciprocal roots without iteration.
    if recover_small_roots:
        monic_quadratic = (
            line_velocity_squared
            + scaled_impact_parameter_squared
            - velocity_intercept_squared
        )
        monic_linear = (
            -2.0 * scaled_line_velocity * scaled_impact_parameter_squared
        )
        monic_constant = line_velocity_squared * scaled_impact_parameter_squared
        if monic_constant != 0.0:
            reciprocal_roots = depressed_quartic(
                monic_constant,
                monic_linear,
                monic_quadratic,
                -2.0 * scaled_line_velocity,
                1.0,
            )
            selected_roots = np.empty(4, dtype=np.float64)
            selected_roots[:] = math.nan
            selected_root_count = 0

            for root in direct_roots:
                if (
                    math.isfinite(root)
                    and abs(root) >= _RECIPROCAL_ROOT_THRESHOLD
                    and selected_root_count < 4
                ):
                    selected_roots[selected_root_count] = root
                    selected_root_count += 1

            for reciprocal_root in reciprocal_roots:
                if math.isfinite(reciprocal_root) and reciprocal_root != 0.0:
                    root = 1.0 / reciprocal_root
                    if (
                        abs(root) < _RECIPROCAL_ROOT_THRESHOLD
                        and selected_root_count < 4
                    ):
                        selected_roots[selected_root_count] = root
                        selected_root_count += 1

            if selected_root_count == direct_real_root_count:
                first_root = selected_roots[0]
                second_root = selected_roots[1]
                third_root = selected_roots[2]
                fourth_root = selected_roots[3]

    return (
        variable_scale * first_root,
        variable_scale * second_root,
        variable_scale * third_root,
        variable_scale * fourth_root,
    )
