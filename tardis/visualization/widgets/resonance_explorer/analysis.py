"""Single-shell geometric and production-solver case analysis."""

import math
from dataclasses import dataclass
from itertools import pairwise

from tardis.transport.montecarlo.configuration.constants import (
    CLOSE_LINE_THRESHOLD,
)

from .geometry import (
    calculate_boundary,
    calculate_geometric_turns,
    calculate_projected_gradient,
    calculate_projected_velocity,
)
from .model import ExplorerInputs
from .solver_adapter import solve_candidates


@dataclass(frozen=True)
class TurnRecord:
    """A projected-gradient zero, including zeros outside the active flight."""

    x: float
    r: float
    s: float
    kind: str
    active: bool
    reason: str
    transport_active: bool
    transport_reason: str


@dataclass(frozen=True)
class CandidateRecord:
    """One raw slot with separate production and equation decisions.

    ``accepted`` and ``event_eligible`` reproduce the production sign,
    distance, shell, and event filters. ``physical`` requires independent
    verification of the original unsquared equation.
    """

    slot_idx: int
    scaled_root: float
    x: float | None
    s: float | None
    r: float | None
    projected_velocity: float | None
    gradient: float | None
    residual: float | None
    accepted: bool
    reason: str
    flags: tuple[str, ...]
    physical: bool
    event_eligible: bool
    unsquared_verified: bool
    verified_event_eligible: bool
    unsquared_tolerance: float | None
    packet_gradient: float | None


@dataclass(frozen=True)
class CaseAnalysis:
    """One immutable analysis snapshot for rendering, UI, and export."""

    inputs: ExplorerInputs
    m: float
    q: float
    r_0: float
    mu_0: float
    outer_chord: tuple[float, float]
    boundary_distance: float
    boundary_kind: str
    x_boundary: float
    active_end: float
    x_end: float
    stop_reason: str
    turn_candidates: tuple[TurnRecord, ...]
    turns: tuple[TurnRecord, ...]
    transport_turns: tuple[TurnRecord, ...]
    monotonic_intervals: tuple[tuple[float, float], ...]
    active_velocity_range: tuple[float, float]
    no_crossing_reason: str | None
    solver_route: str
    raw_roots: tuple[float, float, float, float]
    structured_parameters: tuple[tuple[str, float], ...]
    distance_tolerance: float
    solver_provenance: str
    candidates: tuple[CandidateRecord, ...]
    crossings: tuple[CandidateRecord, ...]
    production_nearest_resonance: CandidateRecord | None
    nearest_resonance: CandidateRecord | None
    analytic_contacts: tuple[TurnRecord, ...]
    ray_family: str
    shape_family: str
    result_family: str
    radial_law: str
    stagnation_radius: float | None
    notes: tuple[str, ...]


def _validate_inputs(inputs: ExplorerInputs) -> None:
    quantities = (
        inputs.r_in,
        inputs.r_out,
        inputs.v_in,
        inputs.v_out,
        inputs.b,
        inputs.x_0,
        inputs.w,
    )
    if not all(math.isfinite(value) for value in quantities):
        raise ValueError("Shell, launch, and velocity inputs must be finite")
    if inputs.r_in <= 0 or inputs.r_out <= inputs.r_in:
        raise ValueError("Shell radii must satisfy 0 < r_in < r_out")
    if inputs.b < 0:
        raise ValueError("Impact parameter b must be nonnegative")
    tolerance = 32 * math.ulp(inputs.r_out)
    if (
        inputs.r_0 < inputs.r_in - tolerance
        or inputs.r_0 > inputs.r_out + tolerance
    ):
        raise ValueError("Launch point must lie within the shell")
    if inputs.s_stop is not None and (
        not math.isfinite(inputs.s_stop) or inputs.s_stop < 0
    ):
        raise ValueError("Earlier stop distance must be finite and nonnegative")


def _classify_turns(
    inputs: ExplorerInputs,
    active_end: float,
    boundary_distance: float,
    distance_tolerance: float,
) -> tuple[TurnRecord, ...]:
    records = []
    # Reproduce the packet interval routine's impact-parameter arithmetic.
    # Its epsilon_d gate is distinct from geometric membership in the ray.
    impact_squared = inputs.r_0**2 * (1.0 - inputs.mu_0**2)
    transport_positions: tuple[float, ...] = ()
    if inputs.m != 0 and inputs.q != 0 and impact_squared > 0:
        turn_radius_cubed = -inputs.q * impact_squared / inputs.m
        if turn_radius_cubed > 0:
            turn_radius = turn_radius_cubed ** (1.0 / 3.0)
            if turn_radius > math.sqrt(impact_squared) + distance_tolerance:
                turn_x = math.sqrt(turn_radius**2 - impact_squared)
                transport_positions = (-turn_x, turn_x)
    for x, kind in calculate_geometric_turns(inputs):
        s = x - inputs.x_0
        r = math.hypot(x, inputs.b)
        if r < inputs.r_in or r > inputs.r_out:
            reason = "outside shell"
        elif s <= 0:
            reason = "behind or at launch"
        elif s >= active_end:
            reason = "at or after active end"
        else:
            reason = "active"
        active = reason == "active"
        if kind != "reversal":
            transport_reason = "stationary inflection does not split transport"
        elif not transport_positions:
            transport_reason = "transport tangent or degenerate tolerance"
        else:
            transport_x = transport_positions[0 if x < 0 else 1]
            transport_s = transport_x - inputs.r_0 * inputs.mu_0
            if transport_s <= distance_tolerance:
                transport_reason = "within launch tolerance"
            elif transport_s >= boundary_distance - distance_tolerance:
                transport_reason = "within boundary tolerance"
            elif transport_s >= active_end:
                transport_reason = "at or after active stop"
            elif not active:
                transport_reason = reason
            else:
                transport_reason = "transport interval split"
        records.append(
            TurnRecord(
                x,
                r,
                s,
                kind,
                active,
                reason,
                transport_reason == "transport interval split",
                transport_reason,
            )
        )
    return tuple(records)


def _classify_candidate(
    slot_idx: int,
    root: float,
    inputs: ExplorerInputs,
    active_end: float,
    boundary_distance: float,
    distance_tolerance: float,
) -> CandidateRecord:
    if not math.isfinite(root):
        return CandidateRecord(
            slot_idx,
            root,
            None,
            None,
            None,
            None,
            None,
            None,
            False,
            "non-real slot",
            ("non-real slot",),
            False,
            False,
            False,
            False,
            None,
            None,
        )
    x = inputs.delta_r * root
    s = x - inputs.x_0
    r = math.hypot(x, inputs.b)
    # The packet filter reconstructs radius from launch state and distance.
    # This can differ by one ulp from hypot(x, b) at a narrow shell edge.
    packet_radius = math.sqrt(
        inputs.r_0 * inputs.r_0 + s * s + 2.0 * inputs.r_0 * s * inputs.mu_0
    )
    flags: list[str] = []
    if r == 0:
        velocity = gradient = packet_gradient = residual = None
        flags.append("coordinate origin undefined")
    else:
        velocity = calculate_projected_velocity(x, inputs)
        gradient = calculate_projected_gradient(x, inputs)
        # Match calculate_packet_position followed by the transport gradient
        # operation order. At tiny G, hypot(x,b) and reconstructed packet
        # position can differ in their last bits.
        packet_mu = (inputs.r_0 * inputs.mu_0 + s) / packet_radius
        packet_radial_velocity = inputs.v_in + inputs.m * (
            packet_radius - inputs.r_in
        )
        packet_gradient = (
            packet_mu**2 * inputs.m
            + (1.0 - packet_mu**2) * packet_radial_velocity / packet_radius
        )
        residual = velocity - inputs.w
    if inputs.delta_v != 0:
        offset = root - inputs.w / inputs.delta_v
        intercept_projection = inputs.q / inputs.delta_v * root
        if (offset > 0 and intercept_projection > 0) or (
            offset < 0 and intercept_projection < 0
        ):
            flags.append("squared-sign artifact")
    # Bound roundoff in m*x + q*x/hypot(x,b) - w by the absolute terms in
    # that unsquared equation. The factor 64 covers the few arithmetic and
    # hypot operations and coefficient rounding, while retaining a bound
    # proportional to local velocity scale instead of an arbitrary fraction
    # of the shell-edge speeds. Near a simple root, dividing this budget by
    # |G| gives the corresponding position uncertainty; at a contact the
    # residual remains the meaningful forward-error criterion.
    unsquared_tolerance = None
    unsquared_verified = False
    if residual is not None:
        radial_term = inputs.q * x / r
        unsquared_tolerance = (
            64
            * math.ulp(1.0)
            * (abs(inputs.m * x) + abs(radial_term) + abs(inputs.w))
        )
        unsquared_verified = abs(residual) <= unsquared_tolerance
        if not unsquared_verified:
            flags.append("unsquared residual too large")
    # Production accepts according to sign and distance, without checking
    # the equation residual. Preserve that decision even for inaccurate slots.
    production_sign_valid = not any(
        flag in flags
        for flag in (
            "coordinate origin undefined",
            "squared-sign artifact",
        )
    )
    physical = production_sign_valid and unsquared_verified
    if s <= distance_tolerance:
        flags.append("at or behind start")
    if s > active_end + distance_tolerance:
        flags.append("after event stop")
    if (
        packet_radius < inputs.r_in - distance_tolerance
        or packet_radius > inputs.r_out + distance_tolerance
    ):
        flags.append("outside shell")
    accepted = production_sign_valid and not any(
        flag in flags
        for flag in (
            "at or behind start",
            "after event stop",
            "outside shell",
        )
    )
    if accepted and s > active_end:
        flags.append("past strict caller end")
    if accepted and s >= boundary_distance and active_end == boundary_distance:
        flags.append("boundary tie")
    event_eligible = accepted and not any(
        flag in flags for flag in ("past strict caller end", "boundary tie")
    )
    verified_event_eligible = event_eligible and physical
    if packet_gradient == 0.0:
        flags.append("singular projected gradient")
    return CandidateRecord(
        slot_idx,
        root,
        x,
        s,
        r,
        velocity,
        gradient,
        residual,
        accepted,
        flags[0] if flags else "eligible resonance",
        tuple(flags),
        physical,
        event_eligible,
        unsquared_verified,
        verified_event_eligible,
        unsquared_tolerance,
        packet_gradient,
    )


def _distinct_crossings(
    candidates: tuple[CandidateRecord, ...],
    turns: tuple[TurnRecord, ...],
    inputs: ExplorerInputs,
    active_end: float,
) -> tuple[CandidateRecord, ...]:
    physical = sorted(
        (
            candidate
            for candidate in candidates
            if candidate.physical
            and candidate.s is not None
            and 0 < candidate.s <= active_end
            and candidate.r is not None
            and inputs.r_in <= candidate.r <= inputs.r_out
        ),
        key=lambda candidate: candidate.x,
    )
    unique: list[CandidateRecord] = []
    for candidate in physical:
        if unique and candidate.x == unique[-1].x:
            continue
        # Merge close slots only when an analytic G=0 contact certifies
        # multiplicity; preserve close simple roots.
        if unique and any(
            turn.active
            and abs(candidate.x - turn.x) <= 1e-7
            and abs(unique[-1].x - turn.x) <= 1e-7
            and calculate_projected_velocity(turn.x, inputs) == inputs.w
            for turn in turns
        ):
            continue
        unique.append(candidate)
    return tuple(unique)


def analyze_case(inputs: ExplorerInputs) -> CaseAnalysis:
    """Analyze a directed shell flight and its production resonance candidates.

    Parameters
    ----------
    inputs : ExplorerInputs
        Shell, ray, target, and optional earlier stop.

    Returns
    -------
    CaseAnalysis
        Immutable geometry, solver, and event-decision snapshot.
    """
    _validate_inputs(inputs)
    boundary_distance, boundary_kind = calculate_boundary(inputs)
    if inputs.s_stop is not None and inputs.s_stop > boundary_distance:
        raise ValueError(
            "Earlier stop must lie within the first shell boundary"
        )
    active_end = boundary_distance if inputs.s_stop is None else inputs.s_stop
    stop_reason = (
        boundary_kind
        if inputs.s_stop is None or inputs.s_stop == boundary_distance
        else inputs.stop_cause
    )
    chord_half = math.sqrt(max(0.0, inputs.r_out**2 - inputs.b**2))
    distance_tolerance = CLOSE_LINE_THRESHOLD * inputs.delta_r
    turn_candidates = _classify_turns(
        inputs, active_end, boundary_distance, distance_tolerance
    )
    turns = tuple(
        turn
        for turn in turn_candidates
        if turn.active and turn.kind == "reversal"
    )
    transport_turns = tuple(
        turn for turn in turn_candidates if turn.transport_active
    )
    points = (
        (inputs.x_0,)
        + tuple(turn.x for turn in transport_turns)
        + (inputs.x_0 + active_end,)
    )
    monotonic_intervals = () if active_end == 0 else tuple(pairwise(points))
    velocity_values = [
        calculate_projected_velocity(inputs.x_0, inputs),
        calculate_projected_velocity(inputs.x_0 + active_end, inputs),
    ]
    velocity_values.extend(
        calculate_projected_velocity(turn.x, inputs) for turn in turns
    )
    active_velocity_range = (min(velocity_values), max(velocity_values))
    solver = solve_candidates(inputs)
    candidates = tuple(
        _classify_candidate(
            slot_idx,
            root,
            inputs,
            active_end,
            boundary_distance,
            distance_tolerance,
        )
        for slot_idx, root in enumerate(solver.roots)
    )
    crossings = _distinct_crossings(
        candidates, turn_candidates, inputs, active_end
    )
    production_eligible = [
        candidate for candidate in candidates if candidate.event_eligible
    ]
    production_nearest = (
        min(production_eligible, key=lambda candidate: candidate.s)
        if production_eligible
        else None
    )
    eligible = [
        candidate
        for candidate in candidates
        if candidate.verified_event_eligible
    ]
    nearest = (
        min(eligible, key=lambda candidate: candidate.s) if eligible else None
    )
    analytic_contacts = tuple(
        turn
        for turn in turn_candidates
        if turn.active
        and calculate_projected_velocity(turn.x, inputs) == inputs.w
    )
    exact_inner_tangent = inputs.x_0 == 0 and inputs.b == inputs.r_in
    if exact_inner_tangent:
        ray_family = "outbound leg"
    elif inputs.x_0 < 0 and boundary_kind == "inner":
        ray_family = "inner-boundary approach"
    elif inputs.x_0 < 0:
        ray_family = "inner-miss flyby"
    elif boundary_kind == "inner":
        ray_family = "inner-boundary approach"
    else:
        ray_family = "outbound leg"
    if active_end == 0 or (inputs.m == 0 and (inputs.q == 0 or inputs.b == 0)):
        shape_family = "flat"
    elif turns:
        shape_family = "one turn" if len(turns) == 1 else "two turns"
    elif any(
        turn.active and turn.kind == "stationary inflection"
        for turn in turn_candidates
    ):
        shape_family = "stationary inflection"
    else:
        shape_family = "monotone"
    if inputs.m == 0 and inputs.q == 0:
        radial_law = "static"
    elif inputs.m == 0:
        radial_law = "uniform radial speed"
    elif inputs.q == 0:
        radial_law = "radius proportional"
    else:
        radial_law = "offset linear"
    stagnation_radius = -inputs.q / inputs.m if inputs.m != 0 else None
    if stagnation_radius is not None and not (
        inputs.r_in <= stagnation_radius <= inputs.r_out
    ):
        stagnation_radius = None
    if solver.non_isolated:
        result_family = "non-isolated match"
    elif analytic_contacts and any(
        turn.kind == "reversal" for turn in analytic_contacts
    ):
        result_family = (
            "extremum touch plus simple crossing"
            if len(crossings) > 1
            else "extremum touch"
        )
    elif analytic_contacts and any(
        turn.kind == "stationary inflection" for turn in analytic_contacts
    ):
        result_family = "flat-center crossing"
    elif not crossings:
        result_family = "no crossing"
    else:
        result_family = {
            1: "one crossing",
            2: "two crossings",
            3: "three crossings",
        }.get(len(crossings), f"{len(crossings)} crossings")
    no_crossing_reason = None
    if result_family == "no crossing":
        outside_range = not (
            active_velocity_range[0] <= inputs.w <= active_velocity_range[1]
        )
        verified_outside = any(
            candidate.physical
            and candidate.s is not None
            and (candidate.s <= 0 or candidate.s > active_end)
            for candidate in candidates
        )
        if outside_range:
            no_crossing_reason = (
                f"Target w={inputs.w:g} is outside active projected-velocity "
                f"range [{active_velocity_range[0]:g}, "
                f"{active_velocity_range[1]:g}]."
            )
        else:
            no_crossing_reason = (
                f"Target w={inputs.w:g} lies within active projected-velocity "
                f"range [{active_velocity_range[0]:g}, "
                f"{active_velocity_range[1]:g}], but no solver slot "
                "supplies a verified crossing."
            )
        if verified_outside:
            no_crossing_reason += (
                " Verified root(s) lie outside active flight or at an "
                "excluded endpoint."
            )
    notes = []
    if exact_inner_tangent:
        notes.append(
            "Exact inner tangent: production boundary code reports immediate "
            "inner hit; geometric outbound continuation has no active "
            "packet interval in this code path."
        )
    if inputs.s_stop is not None and inputs.s_stop < boundary_distance:
        notes.append(f"Active flight truncated by {inputs.stop_cause}.")
    near_zero_differences = [
        candidate
        for candidate in candidates
        if candidate.packet_gradient is not None
        and candidate.packet_gradient != 0.0
        and abs(candidate.gradient) <= 1e-12 * max(abs(inputs.m), 1)
        and candidate.gradient != candidate.packet_gradient
    ]
    if near_zero_differences:
        first = near_zero_differences[0]
        notes.append(
            f"Slot {first.slot_idx} has analytic G={first.gradient:.6g} "
            f"and packet-form G={first.packet_gradient:.6g}; floating-point "
            "rounding leaves the packet value finite, so the exact "
            "singular-depth condition is not met."
        )
    if analytic_contacts:
        notes.append(
            "Analytic G=0 contacts are distinct from production root acceptance."
        )
    unverified_accepted = sum(
        candidate.accepted and not candidate.unsquared_verified
        for candidate in candidates
    )
    if unverified_accepted:
        notes.append(
            f"{unverified_accepted} production-accepted slot(s) failed "
            "the unsquared-equation check and are excluded from verified "
            "crossings and nearest resonance."
        )
    if solver.non_isolated:
        notes.append("Target matches a continuum, not an isolated root.")
    return CaseAnalysis(
        inputs=inputs,
        m=inputs.m,
        q=inputs.q,
        r_0=inputs.r_0,
        mu_0=inputs.mu_0,
        outer_chord=(-chord_half, chord_half),
        boundary_distance=boundary_distance,
        boundary_kind=boundary_kind,
        x_boundary=inputs.x_0 + boundary_distance,
        active_end=active_end,
        x_end=inputs.x_0 + active_end,
        stop_reason=stop_reason,
        turn_candidates=turn_candidates,
        turns=turns,
        transport_turns=transport_turns,
        monotonic_intervals=monotonic_intervals,
        active_velocity_range=active_velocity_range,
        no_crossing_reason=no_crossing_reason,
        solver_route=solver.route,
        raw_roots=solver.roots,
        structured_parameters=solver.structured_parameters,
        distance_tolerance=distance_tolerance,
        solver_provenance=solver.provenance,
        candidates=candidates,
        crossings=crossings,
        production_nearest_resonance=production_nearest,
        nearest_resonance=nearest,
        analytic_contacts=analytic_contacts,
        ray_family=ray_family,
        shape_family=shape_family,
        result_family=result_family,
        radial_law=radial_law,
        stagnation_radius=stagnation_radius,
        notes=tuple(notes),
    )
