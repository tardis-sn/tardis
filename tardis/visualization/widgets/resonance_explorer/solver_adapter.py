"""Production structured quartic and its constant-speed branch."""

import math
from dataclasses import dataclass

from tardis.transport.montecarlo.nonhomologous_grid import (
    solve_resonance_quartic,
)

from .model import ExplorerInputs


@dataclass(frozen=True)
class SolverResult:
    """Four raw production-order root slots and the exact top-level route."""

    route: str
    roots: tuple[float, float, float, float]
    structured_parameters: tuple[tuple[str, float], ...]
    non_isolated: bool
    provenance: str


def solve_candidates(inputs: ExplorerInputs) -> SolverResult:
    """Return structured roots or the production-equivalent analytic limit.

    Returns
    -------
    SolverResult
        Raw dimensionless positions and route metadata.
    """
    delta_r = inputs.delta_r
    delta_v = inputs.delta_v
    q = inputs.q
    w = inputs.w
    nan = math.nan
    transverse_direction_fraction = 1.0 - inputs.mu_0**2
    if delta_v == 0:
        non_isolated = (q == 0 and w == 0) or (
            transverse_direction_fraction == 0
            and w == math.copysign(1.0, inputs.x_0) * q
        )
        if q != 0 and abs(w) < abs(q):
            impact_squared = (
                inputs.r_0 * inputs.r_0 * transverse_direction_fraction
            )
            root_squared = w * w * impact_squared / (q * q - w * w)
            root = math.copysign(math.sqrt(root_squared), w / q) / delta_r
            roots = (root, nan, nan, nan)
        else:
            roots = (nan, nan, nan, nan)
        return SolverResult(
            "Δv=0 analytic",
            roots,
            (),
            non_isolated,
            "code-equivalent constant-speed analytic branch",
        )

    line_velocity = w / delta_v
    scaled_radius = inputs.r_0 / delta_r
    impact_squared = scaled_radius**2 * transverse_direction_fraction
    intercept = q / delta_v
    if impact_squared == 0:
        route = "radial H²=0"
    elif intercept == 0:
        route = "zero-intercept Q=0"
    elif line_velocity == 0:
        route = "line-center N=0"
    else:
        route = "general quartic"
    roots = tuple(
        float(root)
        for root in solve_resonance_quartic(
            line_velocity, impact_squared, intercept
        )
    )
    return SolverResult(
        route,
        roots,
        (("N", line_velocity), ("H²", impact_squared), ("Q", intercept)),
        False,
        "TARDIS solve_resonance_quartic",
    )
