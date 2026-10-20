"""Directed ray and projected-velocity geometry for one shell."""

import math

from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
)

from .model import ExplorerInputs


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
        inputs.r_0, inputs.mu_0, inputs.r_in, inputs.r_out
    )
    return float(distance), "inner" if shell_step == -1 else "outer"


def calculate_projected_velocity(x: float, inputs: ExplorerInputs) -> float:
    """Evaluate the signed packet-direction velocity at ray coordinate ``x``.

    Returns
    -------
    float
        Projected velocity in illustrative velocity units.
    """
    radius = math.hypot(x, inputs.b)
    return inputs.m * x + inputs.q * x / radius


def calculate_projected_gradient(x: float, inputs: ExplorerInputs) -> float:
    """Evaluate the derivative of projected velocity along the ray.

    Returns
    -------
    float
        Velocity per length in illustrative units.
    """
    radius = math.hypot(x, inputs.b)
    return inputs.m + inputs.q * inputs.b**2 / radius**3


def calculate_geometric_turns(
    inputs: ExplorerInputs,
) -> tuple[tuple[float, str], ...]:
    """Find every real geometric zero of the projected gradient.

    Returns
    -------
    tuple of tuple
        Signed ray coordinate and zero type, in travel order.
    """
    if inputs.b == 0 or inputs.m == 0:
        return ()
    if inputs.q == -inputs.m * inputs.b:
        return ((0.0, "stationary inflection"),)
    turn_radius_cubed = -inputs.q * inputs.b**2 / inputs.m
    if turn_radius_cubed <= 0:
        return ()
    turn_radius = turn_radius_cubed ** (1 / 3)
    radius_difference = turn_radius**2 - inputs.b**2
    if radius_difference <= 0:
        return ()
    position = math.sqrt(radius_difference)
    return ((-position, "reversal"), (position, "reversal"))
