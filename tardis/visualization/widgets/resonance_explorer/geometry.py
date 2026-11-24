"""Directed ray and projected-velocity geometry for one shell."""

import math

from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
)
from tardis.visualization.widgets.resonance_explorer.model import ExplorerInputs


def calculate_boundary(inputs: ExplorerInputs) -> tuple[float, str]:
    """Return the production boundary distance and destination.

    Returns
    -------
    distance : float
        Forward distance in illustrative length units.
    boundary : str
        ``inner`` or ``outer`` according to the production decision.
    """
    distance, shell_step = calculate_distance_boundary(
        inputs.r_0, inputs.mu_0, inputs.r_inner, inputs.r_outer
    )
    return float(distance), "inner" if shell_step == -1 else "outer"


def calculate_projected_velocity(x: float, inputs: ExplorerInputs) -> float:
    """Evaluate the signed packet-direction velocity at ray coordinate ``x``.

    Returns
    -------
    float
        Projected velocity in illustrative velocity units.
    """
    radius = math.hypot(x, inputs.impact_parameter)
    return inputs.velocity_slope * x + inputs.velocity_intercept * x / radius


def calculate_projected_gradient(x: float, inputs: ExplorerInputs) -> float:
    """Evaluate the derivative of projected velocity along the ray.

    Returns
    -------
    float
        Velocity per length in illustrative units.
    """
    radius = math.hypot(x, inputs.impact_parameter)
    return (
        inputs.velocity_slope
        + inputs.velocity_intercept * inputs.impact_parameter**2 / radius**3
    )


def calculate_geometric_turns(
    inputs: ExplorerInputs,
) -> tuple[tuple[float, str], ...]:
    """Find every real geometric zero of the projected gradient.

    Returns
    -------
    tuple of tuple
        Signed ray coordinate and zero type, in travel order.
    """
    if inputs.impact_parameter == 0 or inputs.velocity_slope == 0:
        return ()
    if (
        inputs.velocity_intercept
        == -inputs.velocity_slope * inputs.impact_parameter
    ):
        return ((0.0, "stationary inflection"),)
    turn_radius_cubed = (
        -inputs.velocity_intercept
        * inputs.impact_parameter**2
        / inputs.velocity_slope
    )
    if turn_radius_cubed <= 0:
        return ()
    turn_radius = turn_radius_cubed ** (1 / 3)
    radius_difference = turn_radius**2 - inputs.impact_parameter**2
    if radius_difference <= 0:
        return ()
    position = math.sqrt(radius_difference)
    return ((-position, "reversal"), (position, "reversal"))
