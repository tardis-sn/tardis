"""Production structured quartic and its constant-speed branch."""

import math
from dataclasses import dataclass

from tardis.transport.montecarlo.nonhomologous_grid import (
    solve_resonance_quartic,
)
from tardis.visualization.widgets.resonance_explorer.model import ExplorerInputs


@dataclass(frozen=True)
class SolverResult:
    """Four raw production-order root slots and the exact top-level route."""

    route: str
    roots: tuple[float, float, float, float]
    non_isolated: bool


def solve_candidates(inputs: ExplorerInputs) -> SolverResult:
    """Return structured roots or the production-equivalent analytic limit.

    Returns
    -------
    SolverResult
        Raw dimensionless positions and route metadata.
    """
    delta_r = inputs.delta_r
    delta_v = inputs.delta_v
    velocity_intercept = inputs.velocity_intercept
    line_velocity = inputs.line_velocity
    nan = math.nan
    transverse_direction_fraction = 1.0 - inputs.mu_0**2
    if delta_v == 0:
        non_isolated = (velocity_intercept == 0 and line_velocity == 0) or (
            transverse_direction_fraction == 0
            and line_velocity
            == math.copysign(1.0, inputs.launch_position) * velocity_intercept
        )
        if velocity_intercept != 0 and abs(line_velocity) < abs(
            velocity_intercept
        ):
            impact_squared = (
                inputs.r_0 * inputs.r_0 * transverse_direction_fraction
            )
            root_squared = (
                line_velocity
                * line_velocity
                * impact_squared
                / (
                    velocity_intercept * velocity_intercept
                    - line_velocity * line_velocity
                )
            )
            root = (
                math.copysign(
                    math.sqrt(root_squared), line_velocity / velocity_intercept
                )
                / delta_r
            )
            roots = (root, nan, nan, nan)
        else:
            roots = (nan, nan, nan, nan)
        return SolverResult("Δv=0 analytic", roots, non_isolated)

    scaled_line_velocity = line_velocity / delta_v
    scaled_radius = inputs.r_0 / delta_r
    impact_squared = scaled_radius**2 * transverse_direction_fraction
    scaled_velocity_intercept = velocity_intercept / delta_v
    if impact_squared == 0:
        route = "radial H²=0"
    elif scaled_velocity_intercept == 0:
        route = "zero-intercept Q=0"
    elif scaled_line_velocity == 0:
        route = "line-center N=0"
    else:
        route = "general quartic"
    roots = tuple(
        float(root)
        for root in solve_resonance_quartic(
            scaled_line_velocity, impact_squared, scaled_velocity_intercept
        )
    )
    return SolverResult(route, roots, False)
