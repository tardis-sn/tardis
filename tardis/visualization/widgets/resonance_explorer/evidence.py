"""Independent original-equation checks for one explorer snapshot.

The bracket solver here never reads quartic slots or the explorer's interval
classification. It is evidence about the ray equation, not a replacement for
the production transport decision.
"""

import hashlib
import importlib
import math
import subprocess
from dataclasses import asdict
from itertools import pairwise
from pathlib import Path

import numpy as np

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_line_nonhomologous,
)
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    MISS_DISTANCE,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket

from .analysis import CaseAnalysis
from .model import ExplorerInputs

LENGTH_SCALE_CM = 1e14
VELOCITY_SCALE_CM_S = 1e8
PACKET_FREQUENCY_HZ = 1e15


def _source_revision() -> dict[str, object]:
    """Identify the exact explorer and transport source used for evidence."""
    repository_root = Path(__file__).resolve().parents[3]
    solver_path = Path(
        importlib.import_module(
            "tardis.transport.montecarlo.nonhomologous_grid"
        ).__file__
    ).resolve()
    distance_path = Path(
        importlib.import_module(
            "tardis.transport.geometry.calculate_distances"
        ).__file__
    ).resolve()

    def git_commit(root: Path) -> str | None:
        try:
            process = subprocess.run(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                capture_output=True,
                check=False,
                text=True,
                timeout=5,
            )
        except OSError, subprocess.TimeoutExpired:
            return None
        return process.stdout.strip() if process.returncode == 0 else None

    package_dir = Path(__file__).resolve().parent
    return {
        "explorer_commit": git_commit(repository_root),
        "tardis_commit": git_commit(solver_path.parents[3]),
        "explorer_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(package_dir.glob("*.py"))
        },
        "transport_sha256": {
            str(path.relative_to(solver_path.parents[3])): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in (solver_path, distance_path)
        },
    }


def _original_residual(x: float, inputs: ExplorerInputs) -> float:
    """Evaluate the unsquared line equation in illustrative velocity units."""
    radius = math.hypot(x, inputs.b)
    if radius == 0.0:
        return math.nan
    return inputs.m * x + inputs.q * x / radius - inputs.w


def _bracketed_zero(
    lower: float, upper: float, inputs: ExplorerInputs
) -> float:
    """Bisect a strict sign change of the original equation."""
    lower_residual = _original_residual(lower, inputs)
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        if midpoint in (lower, upper):
            break
        midpoint_residual = _original_residual(midpoint, inputs)
        if lower_residual * midpoint_residual <= 0.0:
            upper = midpoint
        else:
            lower = midpoint
            lower_residual = midpoint_residual
    return (lower + upper) / 2.0


def _independent_roots(analysis: CaseAnalysis) -> list[dict[str, object]]:
    """Bracket isolated original-equation roots on independently split ray.

    The derivative zero of ``m*x + q*x/hypot(x,b)`` has radius satisfying
    ``r³=-q*b²/m``. That identity supplies only bracket boundaries; all
    reported roots are checked by direct substitution in the original equation.
    """
    inputs = analysis.inputs
    lower = inputs.x_0
    upper = analysis.x_end
    if upper <= lower or (
        inputs.m == 0
        and _original_residual(lower, inputs) == 0
        and _original_residual(upper, inputs) == 0
    ):
        return []
    points = [lower, upper]
    if inputs.b > 0 and inputs.m != 0 and inputs.m * inputs.q < 0:
        turning_radius_cubed = -inputs.q * inputs.b**2 / inputs.m
        turning_radius = turning_radius_cubed ** (1.0 / 3.0)
        if turning_radius > inputs.b:
            turn_x = math.sqrt(turning_radius**2 - inputs.b**2)
            points.extend(x for x in (-turn_x, turn_x) if lower < x < upper)
    points.sort()
    roots: list[float] = []
    for x in points:
        residual = _original_residual(x, inputs)
        scale = abs(inputs.m * x) + abs(inputs.q) + abs(inputs.w)
        if (
            math.isfinite(residual)
            and abs(residual) <= 32 * math.ulp(1.0) * scale
        ):
            roots.append(x)
    for left, right in pairwise(points):
        left_residual = _original_residual(left, inputs)
        right_residual = _original_residual(right, inputs)
        if left_residual * right_residual < 0.0:
            roots.append(_bracketed_zero(left, right, inputs))
    unique_roots: list[float] = []
    for x in sorted(roots):
        if not unique_roots or abs(x - unique_roots[-1]) > 1e-10:
            unique_roots.append(x)
    result: list[dict[str, object]] = []
    for x in unique_roots:
        s = x - inputs.x_0
        gradient = (
            inputs.m + inputs.q * inputs.b**2 / math.hypot(x, inputs.b) ** 3
        )
        result.append(
            {
                "x": x,
                "s": s,
                "residual": _original_residual(x, inputs),
                "kind": "contact" if abs(gradient) < 1e-10 else "crossing",
                "inside_open_active_flight": 0 < s < analysis.active_end,
            }
        )
    return result


def _packet_nearest(analysis: CaseAnalysis) -> dict[str, object]:
    """Compare with the actual packet distance routine under explicit units."""
    inputs = analysis.inputs
    velocities = (inputs.v_in, inputs.v_out, inputs.w)
    if max(abs(value) for value in velocities) > 10 or inputs.r_out > 1e4:
        return {
            "status": "not mapped",
            "reason": "Illustrative input exceeds the documented physical-unit mapping range.",
        }
    if analysis.result_family == "non-isolated match":
        return {
            "status": "not mapped",
            "reason": "A continuum has no isolated nearest line distance.",
        }
    geometry = NumbaRadial1DGeometry(
        np.asarray([inputs.r_in * LENGTH_SCALE_CM]),
        np.asarray([inputs.r_out * LENGTH_SCALE_CM]),
        np.asarray([inputs.v_in * VELOCITY_SCALE_CM_S]),
        np.asarray([inputs.v_out * VELOCITY_SCALE_CM_S]),
    )
    packet = RPacket(
        inputs.r_0 * LENGTH_SCALE_CM,
        inputs.mu_0,
        PACKET_FREQUENCY_HZ,
        1.0,
        1963,
    )
    line_frequency = float(
        PACKET_FREQUENCY_HZ
        * (1.0 - inputs.w * VELOCITY_SCALE_CM_S / C_SPEED_OF_LIGHT)
    )
    recovered_target = float(
        C_SPEED_OF_LIGHT
        * (1.0 - line_frequency / PACKET_FREQUENCY_HZ)
        / VELOCITY_SCALE_CM_S
    )
    if abs(recovered_target - inputs.w) > 1e-10:
        return {
            "status": "not mapped",
            "reason": "Frequency rounding cannot represent the target within 1e-10 velocity units.",
        }
    packet_distance_cm = calculate_distance_line_nonhomologous(
        packet,
        geometry,
        line_frequency,
        maximum_distance=analysis.active_end * LENGTH_SCALE_CM,
    )
    packet_distance = (
        None
        if packet_distance_cm == MISS_DISTANCE
        else packet_distance_cm / LENGTH_SCALE_CM
    )
    explorer_distance = (
        None
        if analysis.production_nearest_resonance is None
        else analysis.production_nearest_resonance.s
    )
    difference = (
        None
        if packet_distance is None or explorer_distance is None
        else packet_distance - explorer_distance
    )
    contact_shifted = bool(analysis.analytic_contacts) and (
        recovered_target != inputs.w
    )
    return {
        "status": (
            "contact sensitive to frequency rounding"
            if contact_shifted
            else "compared"
            if packet_distance is not None or explorer_distance is not None
            else "both miss"
        ),
        "comparison_valid": not contact_shifted,
        "reason": (
            "The physical-frequency mapping shifts an exact tangent target; "
            "the resulting packet distance is displayed but is not exact "
            "parity evidence for the illustrative contact."
            if contact_shifted
            else None
        ),
        "length_scale_cm": LENGTH_SCALE_CM,
        "velocity_scale_cm_s": VELOCITY_SCALE_CM_S,
        "packet_frequency_hz": PACKET_FREQUENCY_HZ,
        "line_frequency_hz": line_frequency,
        "recovered_target": recovered_target,
        "packet_distance": packet_distance,
        "explorer_production_distance": explorer_distance,
        "difference": difference,
        "agreement_within_1e-8": (
            (packet_distance is None and explorer_distance is None)
            or (difference is not None and abs(difference) <= 1e-8)
        ),
    }


def build_evidence_record(analysis: CaseAnalysis) -> dict[str, object]:
    """Build independent equation and packet-distance evidence for a case.

    Parameters
    ----------
    analysis : CaseAnalysis
        Current immutable explorer snapshot.

    Returns
    -------
    dict[str, object]
        JSON-safe original-equation roots, raw production slots, and mapped
        packet distance. Source revision is carried by export provenance.
    """
    independent_roots = _independent_roots(analysis)
    production_slots = [
        {
            "slot_idx": candidate.slot_idx,
            "scaled_root": (
                candidate.scaled_root
                if math.isfinite(candidate.scaled_root)
                else None
            ),
            "disposition": candidate.reason,
            "flags": list(candidate.flags),
            "production_accepted": candidate.accepted,
            "unsquared_verified": candidate.unsquared_verified,
            "event_eligible": candidate.event_eligible,
            "residual": candidate.residual,
        }
        for candidate in analysis.candidates
    ]
    unmatched = [
        root["x"]
        for root in independent_roots
        if root["inside_open_active_flight"]
        and not any(
            candidate.x is not None
            and candidate.unsquared_verified
            and abs(candidate.x - root["x"]) <= 1e-7
            for candidate in analysis.candidates
        )
    ]
    return {
        "method": "original unsquared equation, derivative-bounded bisection",
        "inputs": asdict(analysis.inputs),
        "route": analysis.solver_route,
        "independent_roots": independent_roots,
        "production_slots": production_slots,
        "agreement": {
            "status": "mismatch" if unmatched else "matched",
            "unmatched_independent_roots": unmatched,
        },
        "packet_nearest": _packet_nearest(analysis),
        "source_revision": _source_revision(),
    }
