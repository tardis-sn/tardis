"""Linked Bokeh views of one analyzed single-shell packet flight.

The views are built once and then refilled from each analysis snapshot, so a
slider change sends only new plot data to the browser. The renderer samples
the displayed equations but does not classify roots, turns, or event
ownership.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl
import numpy as np
import numpy.typing as npt
import panel as pn
from bokeh.layouts import column, row
from bokeh.models import (
    Arrow,
    ColumnDataSource,
    Label,
    Legend,
    LegendItem,
    Model,
    NormalHead,
    Range1d,
)
from bokeh.plotting import figure

if TYPE_CHECKING:
    from bokeh.models import LayoutDOM

    from tardis.visualization.widgets.resonance_explorer.analysis import (
        CaseAnalysis,
    )

INK = "#30343b"
MUTED = "#899198"
RED = "#b44943"
GREEN = "#3b7667"
AMBER = "#a66b2e"
UNVERIFIED = "#a44d77"
STOP = "#74487c"
SHELL_PALETTE = [
    mpl.colors.to_hex(color)
    for color in mpl.colormaps["RdBu_r"](np.linspace(0.20, 0.80, 256))
]
DASHED = "4 3"
SHORT_DASHED = "3 3"
DOTTED = "1 4"

# Plot frame sizes in screen pixels
PROJECTED_WIDTH = 620
PROJECTED_HEIGHT = 400
GRADIENT_HEIGHT = 110
SIDE_WIDTH = 220
RADIAL_HEIGHT = 110
# Side panels share one left border so Bokeh aligns their frames without
# pushing the axis-free panels past their own width
SIDE_BORDER = 50

# Screen pixels per typographic point, to keep label offsets readable
PIXELS_PER_POINT = 4 / 3

MARKER_COLUMNS = (
    "x",
    "y",
    "marker",
    "size",
    "angle",
    "fill_color",
    "fill_alpha",
    "line_color",
    "line_width",
    "label",
)
LABEL_COLUMNS = (
    "x",
    "y",
    "text",
    "x_offset",
    "y_offset",
    "text_color",
    "text_font_size",
    "text_align",
    "text_baseline",
)
GUIDE_COLUMNS = (
    "x",
    "line_color",
    "line_dash",
    "line_alpha",
    "line_width",
    "label",
)
SEGMENT_COLUMNS = (
    "x0",
    "y0",
    "x1",
    "y1",
    "line_color",
    "line_width",
    "line_dash",
    "line_alpha",
)
ARROW_COLUMNS = (
    "x_start",
    "y_start",
    "x_end",
    "y_end",
    "line_color",
    "line_width",
    "head_size",
    "head_color",
    "head_line_color",
)


def _columns(rows: list[dict], keys: tuple[str, ...]) -> dict[str, list]:
    """Turn plot rows into source columns, keeping every column when empty."""
    return {key: [item[key] for item in rows] for key in keys}


def _marker(
    x: float,
    y: float,
    marker: str,
    size: float,
    line_color: str,
    *,
    fill_color: str | None = "white",
    line_width: float = 1.4,
    angle: float = 0.0,
    label: str = "",
    z: float = 0.0,
) -> dict:
    """Describe one scatter marker; a fill colour of None leaves it hollow."""
    return dict(
        x=x,
        y=y,
        marker=marker,
        size=size,
        angle=angle,
        fill_color=fill_color or "white",
        fill_alpha=0.0 if fill_color is None else 1.0,
        line_color=line_color,
        line_width=line_width,
        label=label,
        z=z,
    )


def _label(
    x: float,
    y: float,
    text: str,
    *,
    offset: tuple[float, float] = (0.0, 0.0),
    color: str = INK,
    size: float = 8.0,
    align: str = "center",
    baseline: str = "bottom",
) -> dict:
    """Describe one text label offset in points, upward positive."""
    return dict(
        x=x,
        y=y,
        text=text,
        x_offset=offset[0] * PIXELS_PER_POINT,
        y_offset=-offset[1] * PIXELS_PER_POINT,
        text_color=color,
        text_font_size=f"{size}pt",
        text_align=align,
        text_baseline=baseline,
    )


def _segment(
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    width: float,
    *,
    dash: str = "solid",
    alpha: float = 1.0,
) -> dict:
    """Describe one straight line segment."""
    return dict(
        x0=start[0],
        y0=start[1],
        x1=end[0],
        y1=end[1],
        line_color=color,
        line_width=width,
        line_dash=dash,
        line_alpha=alpha,
    )


def _arrow(
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    width: float,
    *,
    head_size: float = 8.0,
    head_line_color: str | None = None,
) -> dict:
    """Describe one arrow with a filled head."""
    return dict(
        x_start=start[0],
        y_start=start[1],
        x_end=end[0],
        y_end=end[1],
        line_color=color,
        line_width=width,
        head_size=head_size,
        head_color=color,
        head_line_color=head_line_color or color,
    )


def _shell_color(speed: float, speed_scale: float) -> str:
    """Map a signed radial speed onto the symmetric shell palette."""
    fraction = min(max(0.5 + 0.5 * speed / speed_scale, 0.0), 1.0)
    return SHELL_PALETTE[min(int(fraction * len(SHELL_PALETTE)), 255)]


def _sample_flight(
    coordinates: npt.NDArray[np.float64], analysis: CaseAnalysis
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Sample v_parallel and G on the ray, blank outside the shell."""
    impact = analysis.inputs.impact_parameter
    radius = np.hypot(coordinates, impact)
    projected = np.full_like(coordinates, np.nan)
    gradient = np.full_like(coordinates, np.nan)
    np.divide(
        analysis.velocity_intercept * coordinates,
        radius,
        out=projected,
        where=radius > 0,
    )
    projected += analysis.velocity_slope * coordinates
    np.divide(
        analysis.velocity_intercept * impact**2,
        radius**3,
        out=gradient,
        where=radius > 0,
    )
    gradient += analysis.velocity_slope
    in_shell = (radius >= analysis.inputs.r_inner) & (
        radius <= analysis.inputs.r_outer
    )
    projected[~in_shell] = np.nan
    gradient[~in_shell] = np.nan
    return projected, gradient


def _padded_limits(values: list[float]) -> tuple[float, float]:
    """Pad a value span so markers on its edges stay inside the frame."""
    lower, upper = float(np.min(values)), float(np.max(values))
    padding = max(0.08 * (upper - lower), 0.12)
    return lower - padding, upper + padding


def _turn_color(turn: object) -> str:
    """Colour a gradient zero by whether transport passes through it."""
    if turn.transport_active or (turn.active and "inflection" in turn.kind):
        return GREEN
    return AMBER if turn.active else MUTED


def _candidate_style(candidate: object) -> tuple[str, str, str]:
    """Choose marker shape, edge, and label from a candidate disposition."""
    if candidate.verified_event_eligible:
        return "circle", RED, "resonance on path"
    if candidate.accepted and not candidate.unsquared_verified:
        return "square_x", UNVERIFIED, "root failing the original equation"
    if candidate.physical:
        return "x", AMBER, "root not reached"
    if candidate.x is not None and np.isfinite(candidate.x):
        return "triangle", MUTED, "wrong-sign root"
    return "", MUTED, "unavailable slot"


def _guide_rows(analysis: CaseAnalysis) -> tuple[list[dict], list[dict]]:
    """Mark launch, edge, and end (legend) and turns and crossings (no legend)."""
    flight_ends_early = analysis.x_end != analysis.x_boundary
    labelled = [
        dict(
            x=x,
            line_color=color,
            line_dash=dash,
            line_alpha=0.65,
            line_width=0.85,
            label=label,
        )
        for x, color, dash, label in (
            (analysis.inputs.launch_position, INK, DASHED, "launch position"),
            (analysis.x_boundary, MUTED, DASHED, "shell edge reached"),
            (
                analysis.x_end,
                STOP if flight_ends_early else MUTED,
                SHORT_DASHED,
                "flight end" if flight_ends_early else "shell edge reached",
            ),
        )
    ]
    unlabelled = [
        dict(
            x=turn.x,
            line_color=_turn_color(turn),
            line_dash=DOTTED,
            line_alpha=0.35,
            line_width=0.65,
            label="",
        )
        for turn in analysis.turn_candidates
        if np.isfinite(turn.x)
    ] + [
        dict(
            x=crossing.x,
            line_color=RED,
            line_dash=DOTTED,
            line_alpha=0.3,
            line_width=0.75,
            label="",
        )
        for crossing in analysis.crossings
        if np.isfinite(crossing.x)
    ]
    return labelled, unlabelled


def _event_rows(
    analysis: CaseAnalysis,
    selected_slot_index: int | None,
    raw_slot_height: float,
    geometry_pixels_per_unit: float,
) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    """Place turns, raw slots, candidates, and numbered crossings per panel."""
    inputs = analysis.inputs
    impact = inputs.impact_parameter
    markers = {"geometry": [], "projected": [], "gradient": []}
    labels = {"geometry": [], "projected": []}
    for turn in analysis.turn_candidates:
        if not np.isfinite(turn.x):
            continue
        inflection = "inflection" in turn.kind
        color = _turn_color(turn)
        fill = color if turn.transport_active and not inflection else "white"
        label = "G = 0 inflection" if inflection else "G = 0"
        for panel, ordinate in (
            ("geometry", impact),
            (
                "projected",
                analysis.velocity_slope * turn.x
                + analysis.velocity_intercept * turn.x / turn.r,
            ),
            ("gradient", 0.0),
        ):
            markers[panel].append(
                _marker(
                    turn.x,
                    ordinate,
                    "square" if inflection else "diamond",
                    9.5,
                    color,
                    fill_color=fill,
                    label=label,
                    z=7,
                )
            )
        if turn in analysis.analytic_contacts:
            for panel, ordinate in (
                ("geometry", impact),
                ("projected", inputs.line_velocity),
            ):
                markers[panel].append(
                    _marker(
                        turn.x,
                        ordinate,
                        "plus" if inflection else "star",
                        15.5,
                        GREEN,
                        line_width=1.6,
                        label="target touches curve",
                        z=9,
                    )
                )

    crossing_coordinates = {crossing.x for crossing in analysis.crossings}
    for candidate in analysis.candidates:
        if candidate.x is None or not np.isfinite(candidate.x):
            continue
        markers["projected"].append(
            _marker(
                candidate.x,
                raw_slot_height,
                "inverted_triangle",
                7.5,
                MUTED,
                line_width=1.1,
                label="quartic root position",
                z=7,
            )
        )
        marker, color, label = _candidate_style(candidate)
        if not marker:
            continue
        # Numbered crossings own the prominent marker; raw-slot artifacts and
        # excluded roots keep their own shapes
        if candidate.x not in crossing_coordinates or not candidate.physical:
            for panel, ordinate in (
                ("geometry", impact),
                ("projected", candidate.projected_velocity),
            ):
                if ordinate is not None and np.isfinite(ordinate):
                    markers[panel].append(
                        _marker(
                            candidate.x,
                            ordinate,
                            marker,
                            8.7,
                            color,
                            line_width=1.2,
                            label=label,
                            z=8,
                        )
                    )
        if candidate.slot_idx == selected_slot_index:
            for panel, ordinate in (
                ("geometry", impact),
                ("projected", candidate.projected_velocity),
                ("gradient", candidate.gradient),
            ):
                if ordinate is not None and np.isfinite(ordinate):
                    markers[panel].append(
                        _marker(
                            candidate.x,
                            ordinate,
                            "circle",
                            18.0,
                            STOP,
                            fill_color=None,
                            line_width=1.7,
                            z=10,
                        )
                    )

    label_positions: list[tuple[float, int]] = []
    for number, crossing in enumerate(analysis.crossings, start=1):
        if crossing.x is None or not np.isfinite(crossing.x):
            continue
        pixel_x = crossing.x * geometry_pixels_per_unit
        label_tier = 0
        while any(
            tier == label_tier and abs(pixel_x - previous_x) < 20
            for previous_x, tier in label_positions
        ):
            label_tier += 1
        label_positions.append((pixel_x, label_tier))
        face = RED if crossing.verified_event_eligible else "white"
        for panel, ordinate in (
            ("geometry", impact),
            ("projected", inputs.line_velocity),
        ):
            markers[panel].append(
                _marker(
                    crossing.x,
                    ordinate,
                    "circle",
                    10.8,
                    RED,
                    fill_color=face,
                    label="resonance on path (numbered)",
                    z=9,
                )
            )
            labels[panel].append(
                _label(
                    crossing.x,
                    ordinate,
                    str(number),
                    offset=(0, 9 + 13 * label_tier),
                    color=RED,
                )
            )
    nearest = analysis.nearest_resonance
    if nearest is not None and nearest.x is not None:
        for panel, ordinate in (
            ("geometry", impact),
            ("projected", inputs.line_velocity),
        ):
            markers[panel].append(
                _marker(
                    nearest.x,
                    ordinate,
                    "circle",
                    16.0,
                    INK,
                    fill_color=None,
                    line_width=1.7,
                    label="first resonance on path",
                    z=9.5,
                )
            )
    for rows in markers.values():
        rows.sort(key=lambda item: item["z"])
    return markers, labels


def _vector_construction(
    analysis: CaseAnalysis,
    limits: tuple[float, float],
    speed_scale: float,
) -> tuple[dict[str, dict], dict[str, dict]]:
    """Draw radial velocities and their signed projections on the packet ray.

    Only distinct physical forward crossings get vectors; algebraic slots do
    not. The panel shares the shell plot limits, so moving the launch point
    cannot shift the plotted origin.
    """
    inputs = analysis.inputs
    impact = inputs.impact_parameter
    chord_start, chord_end = analysis.outer_chord
    chord_span = max(chord_end - chord_start, 0.5)
    vector_data: list[tuple[float, float, float, float]] = []
    for crossing in analysis.crossings:
        if crossing.x is None:
            continue
        radius = np.hypot(crossing.x, impact)
        if radius <= 0:
            continue
        radial_speed = (
            analysis.velocity_slope * radius + analysis.velocity_intercept
        )
        if abs(radial_speed) <= 1e-13 * max(
            1.0,
            abs(analysis.velocity_slope * radius),
            abs(analysis.velocity_intercept),
        ):
            radial_speed = 0.0
        vector_data.append(
            (
                crossing.x,
                radial_speed * crossing.x / radius,
                radial_speed * impact / radius,
                radial_speed,
            )
        )
    largest_speed = max((abs(item[3]) for item in vector_data), default=0.0)
    largest_vertical = max((abs(item[2]) for item in vector_data), default=0.0)
    horizontal_scale = (
        0.32 * chord_span / largest_speed if largest_speed > 0 else np.inf
    )
    vertical_scale = (
        1.28 * max(impact, 0.6) / largest_vertical
        if largest_vertical > 0
        else np.inf
    )
    vector_scale = min(horizontal_scale, vertical_scale)
    if not np.isfinite(vector_scale):
        vector_scale = 1.0
    if impact == 0:
        for x, projected_speed, _, _ in vector_data:
            if projected_speed < 0 and x > inputs.launch_position:
                vector_scale = min(
                    vector_scale,
                    0.7 * (x - inputs.launch_position) / -projected_speed,
                )
            elif projected_speed > 0 and x < analysis.x_end:
                vector_scale = min(
                    vector_scale, 0.7 * (analysis.x_end - x) / projected_speed
                )
    margin = 0.03 * (limits[1] - limits[0])
    for x, projected_speed, vertical_speed, _ in vector_data:
        for position, speed in ((x, projected_speed), (impact, vertical_speed)):
            if speed > 0:
                vector_scale = min(
                    vector_scale, (limits[1] - margin - position) / speed
                )
            elif speed < 0:
                vector_scale = min(
                    vector_scale, (position - limits[0] - margin) / -speed
                )

    start = (inputs.launch_position, impact)
    end = (analysis.x_end, impact)
    segments = [
        _segment((chord_start, impact), (chord_end, impact), MUTED, 1.0),
        _segment(start, end, "#727a81", 1.1),
    ]
    arrows = []
    if analysis.x_end > inputs.launch_position:
        flight_span = analysis.x_end - inputs.launch_position
        arrows.append(
            _arrow(
                (inputs.launch_position + 0.70 * flight_span, impact),
                (inputs.launch_position + 0.76 * flight_span, impact),
                "#727a81",
                1.1,
                head_size=10.0,
            )
        )
    markers = [
        _marker(0.0, 0.0, "circle", 5.0, MUTED, fill_color=MUTED),
        *(
            _marker(
                *point,
                "dash",
                12.0,
                MUTED,
                line_width=1.0,
                angle=np.pi / 2,
            )
            for point in (start, end)
        ),
    ]
    labels = [
        _label(0.0, 0.0, "O", offset=(5, -12), size=9, align="left"),
        _label(*start, "start", offset=(0, 8), color=MUTED),
        _label(*end, "end", offset=(0, -16), color=MUTED),
    ]
    close_numbers: list[tuple[float, int]] = []
    clustered = False
    radial_width = 3.0 if impact == 0 else 1.8
    for number, (x, projected_speed, vertical_speed, radial_speed) in enumerate(
        vector_data, start=1
    ):
        radial_color = _shell_color(radial_speed, speed_scale)
        point = (x, impact)
        radial_tip = (
            x + vector_scale * projected_speed,
            impact + vector_scale * vertical_speed,
        )
        projection_tip = (x + vector_scale * projected_speed, impact)
        segments.append(
            _segment((0.0, 0.0), point, MUTED, 0.9, dash="2 3", alpha=0.55)
        )
        if radial_speed == 0.0:
            markers.append(
                _marker(
                    *point,
                    "circle",
                    12.3,
                    "black",
                    fill_color=radial_color,
                    line_width=0.7,
                )
            )
        else:
            # A black underlay outlines the shell-coloured shaft
            segments.append(
                _segment(point, radial_tip, "black", radial_width + 0.7)
            )
            arrows.append(
                _arrow(
                    point,
                    radial_tip,
                    radial_color,
                    radial_width,
                    head_line_color="black",
                )
            )
        if np.isclose(projected_speed, 0.0, rtol=0, atol=1e-14):
            markers.append(_marker(*point, "circle", 12.6, RED, line_width=1.8))
        else:
            arrows.append(_arrow(point, projection_tip, RED, 1.7))
        segments.append(
            _segment(radial_tip, projection_tip, MUTED, 1.0, dash=DASHED)
        )
        markers.append(
            _marker(
                *point, "circle", 8.8, "white", fill_color=RED, line_width=0.8
            )
        )
        label_tier = 0
        while any(
            tier == label_tier and abs(x - previous_x) < 0.07 * chord_span
            for previous_x, tier in close_numbers
        ):
            label_tier += 1
            clustered = True
        close_numbers.append((x, label_tier))
        labels.append(
            _label(
                *point,
                str(number),
                offset=(0, 11 + 17 * label_tier),
                color=RED,
                size=8.5,
            )
        )
        if radial_speed == 0.0:
            labels.append(
                _label(*point, "vᵣ=0", offset=(10, -14), size=7.5, align="left")
            )
        if projected_speed == 0.0:
            labels.append(
                _label(
                    *point,
                    "v∥=0",
                    offset=(10, 8),
                    color=RED,
                    size=7.5,
                    align="left",
                )
            )

    if impact == 0 and vector_data:
        note = "vᵣ = v∥; arrows overlap"
    elif clustered:
        note = "Nearby resonance labels"
    else:
        note = ""
    sources = {
        "vector-segments": _columns(segments, SEGMENT_COLUMNS),
        "vector-arrows": _columns(arrows, ARROW_COLUMNS),
        "vector-markers": _columns(markers, MARKER_COLUMNS),
        "vector-labels": _columns(labels, LABEL_COLUMNS),
    }
    properties = {
        "vector-empty-note": dict(visible=not vector_data),
        "vector-note": dict(text=note),
        "vector-legend-radial": dict(
            line_color=_shell_color(
                vector_data[0][3] if vector_data else 0.0, speed_scale
            )
        ),
    }
    return sources, properties


def calculate_guided_plot_data(
    primary: CaseAnalysis,
    reference: CaseAnalysis | None = None,
    selected_slot_index: int | None = None,
    *,
    show_reference: bool = True,
) -> tuple[dict[str, dict], dict[str, dict]]:
    """Compute every plotted value of the guided view.

    Parameters
    ----------
    primary : CaseAnalysis
        Current shell flight and the only source of candidate decisions.
    reference : CaseAnalysis or None, optional
        Paired comparison. It contributes to shared scales even when hidden,
        but never contributes candidate or event markers.
    selected_slot_index : int or None, optional
        Production slot to ring in all applicable views.
    show_reference : bool, optional
        Draw the paired comparison curves when a reference is supplied.

    Returns
    -------
    sources : dict
        Column data for each named ColumnDataSource.
    properties : dict
        Property values for each named range, label, and glyph.
    """
    inputs = primary.inputs
    impact = inputs.impact_parameter
    cases = (primary,) if reference is None else (primary, reference)
    outer_edge = max(
        np.sqrt(
            max(case.inputs.r_outer**2 - case.inputs.impact_parameter**2, 0.0)
        )
        for case in cases
    )
    half_width = max(1.08 * max(case.inputs.r_outer for case in cases), 0.5)
    coordinates = np.linspace(-outer_edge, outer_edge, 1201)
    sampled = [_sample_flight(coordinates, case) for case in cases]
    radii = [
        np.linspace(case.inputs.r_inner, case.inputs.r_outer, 220)
        for case in cases
    ]
    radial_values = [
        case.velocity_slope * case_radii + case.velocity_intercept
        for case, case_radii in zip(cases, radii, strict=True)
    ]
    projected_extrema = [inputs.line_velocity, 0.0]
    gradient_extrema = [0.0]
    radial_extrema = [0.0]
    for (projected, gradient), radial in zip(
        sampled, radial_values, strict=True
    ):
        projected_extrema.extend(projected[np.isfinite(projected)])
        gradient_extrema.extend(gradient[np.isfinite(gradient)])
        radial_extrema.extend(radial)
    projected_limits = _padded_limits(projected_extrema)
    speed_scale = max(abs(value) for value in radial_extrema) or 1e-9

    shown_reference = reference is not None and show_reference
    current_label = f"current {'homology' if primary.velocity_intercept == 0 else 'nonhomology'}"
    reference_label = (
        f"{'homology' if reference.velocity_intercept == 0 else 'nonhomology'} reference"
        if shown_reference
        else ""
    )
    active_coordinates = np.linspace(inputs.launch_position, primary.x_end, 500)
    active_projection, active_gradient = _sample_flight(
        active_coordinates, primary
    )
    reference_rows = slice(None) if shown_reference else slice(0)

    guides, guide_marks = _guide_rows(primary)
    markers, labels = _event_rows(
        primary,
        selected_slot_index,
        raw_slot_height=projected_limits[0]
        + 0.035 * (projected_limits[1] - projected_limits[0]),
        geometry_pixels_per_unit=SIDE_WIDTH / (2 * half_width),
    )
    projected_markers = [row for row in markers["projected"] if row["label"]]
    projected_rings = [row for row in markers["projected"] if not row["label"]]

    ring_radii = np.linspace(inputs.r_outer, inputs.r_inner, 96)
    stagnation = primary.stagnation_radius
    has_stagnation = (
        stagnation is not None
        and inputs.r_inner <= stagnation <= inputs.r_outer
    )
    geometry_labels = [
        _label(*start_or_end, text, offset=(0, 8), size=7)
        for start_or_end, text in (
            ((inputs.launch_position, impact), "start"),
            ((primary.x_end, impact), "end"),
        )
    ] + labels["geometry"]
    if has_stagnation:
        geometry_labels.append(
            _label(
                0.0,
                -stagnation,
                "  radial v=0",
                align="left",
                baseline="middle",
            )
        )

    sources = {
        "active-shading": dict(x0=[inputs.launch_position], x1=[primary.x_end]),
        "line-target": dict(y=[inputs.line_velocity]),
        "projected-curve": dict(x=coordinates, y=sampled[0][0]),
        "projected-active": dict(
            x=active_coordinates,
            y=active_projection,
            label=[current_label] * len(active_coordinates),
        ),
        "projected-reference": dict(
            x=coordinates[reference_rows],
            y=sampled[-1][0][reference_rows],
            label=[reference_label] * len(coordinates[reference_rows]),
        ),
        "gradient-curve": dict(x=coordinates, y=sampled[0][1]),
        "gradient-active": dict(x=active_coordinates, y=active_gradient),
        "gradient-reference": dict(
            x=coordinates[reference_rows], y=sampled[-1][1][reference_rows]
        ),
        "guides": _columns(guides, GUIDE_COLUMNS),
        "guide-marks": _columns(guide_marks, GUIDE_COLUMNS),
        "projected-markers": _columns(projected_markers, MARKER_COLUMNS),
        "projected-rings": _columns(projected_rings, MARKER_COLUMNS),
        "gradient-markers": _columns(markers["gradient"], MARKER_COLUMNS),
        "geometry-markers": _columns(
            [
                _marker(
                    inputs.launch_position,
                    impact,
                    "circle",
                    7.8,
                    INK,
                    fill_color=INK,
                ),
                _marker(primary.x_end, impact, "square", 8.0, INK),
                *markers["geometry"],
            ],
            MARKER_COLUMNS,
        ),
        "projected-labels": _columns(labels["projected"], LABEL_COLUMNS),
        "geometry-labels": _columns(geometry_labels, LABEL_COLUMNS),
        "shell-rings": dict(
            inner_radius=ring_radii[1:],
            outer_radius=ring_radii[:-1],
            fill_color=[
                _shell_color(
                    primary.velocity_slope * radius
                    + primary.velocity_intercept,
                    speed_scale,
                )
                for radius in ring_radii[:-1]
            ],
        ),
        "shell-edges": dict(
            radius=[inputs.r_inner, inputs.r_outer], fill_alpha=[1.0, 0.0]
        ),
        "stagnation-ring": dict(radius=[stagnation] if has_stagnation else []),
        "geometry-ray": dict(y=[impact]),
        "geometry-active-ray": dict(
            x0=[inputs.launch_position], x1=[primary.x_end], y=[impact]
        ),
        "radial-current": dict(
            x=radii[0],
            y=radial_values[0],
            label=[current_label] * len(radii[0]),
        ),
        "radial-reference": dict(
            x=radii[-1][reference_rows],
            y=radial_values[-1][reference_rows],
            label=[reference_label] * len(radii[-1][reference_rows]),
        ),
        "radial-edges": dict(
            x=[inputs.r_inner, inputs.r_outer],
            y=[inputs.v_inner, inputs.v_outer],
        ),
    }
    properties = {
        "projected-x": dict(start=-half_width, end=half_width),
        "projected-y": dict(
            zip(("start", "end"), projected_limits, strict=True)
        ),
        "gradient-y": dict(
            zip(("start", "end"), _padded_limits(gradient_extrema), strict=True)
        ),
        "geometry-x": dict(start=-half_width, end=half_width),
        "geometry-y": dict(start=-half_width, end=half_width),
        "radial-x": dict(
            start=inputs.r_inner - 0.08, end=inputs.r_outer + 0.08
        ),
        "radial-y": dict(
            zip(("start", "end"), _padded_limits(radial_extrema), strict=True)
        ),
        "line-target-label": dict(
            text=f"line target w = {inputs.line_velocity:.6g}"
        ),
        "radial-edges-label": dict(
            text=f"edges: {inputs.v_inner:.5g}, {inputs.v_outer:.5g}"
        ),
    }
    vector_sources, vector_properties = _vector_construction(
        primary, (-half_width, half_width), speed_scale
    )
    # Plain lists travel as JSON; NumPy arrays travel as binary buffers, and
    # buffer updates did not reach the VS Code notebook view of the plot
    sources = {
        name: {
            column: values.tolist()
            if isinstance(values, np.ndarray)
            else values
            for column, values in data.items()
        }
        for name, data in (sources | vector_sources).items()
    }
    return sources, properties | vector_properties


def _figure(
    x_range: Range1d,
    y_range: Range1d,
    width: int,
    height: int,
    *,
    axes: bool = True,
    **kwargs: object,
) -> figure:
    """Create one fixed-size, tool-free panel in the shared style."""
    plot = figure(
        x_range=x_range,
        y_range=y_range,
        frame_width=width,
        frame_height=height,
        tools="",
        toolbar_location=None,
        background_fill_color="white",
        border_fill_color="white",
        outline_line_color=None,
        x_axis_location="below" if axes else None,
        y_axis_location="left" if axes else None,
        **kwargs,
    )
    plot.grid.visible = False
    plot.title.text_font_size = "10pt"
    plot.title.text_font_style = "normal"
    plot.title.text_color = INK
    for axis in plot.axis:
        axis.axis_line_color = MUTED
        axis.major_tick_line_color = MUTED
        axis.minor_tick_line_color = None
        axis.major_label_text_color = INK
        axis.major_label_text_font_size = "8pt"
        axis.axis_label_text_font_size = "9pt"
        axis.axis_label_text_font_style = "normal"
        axis.axis_label_text_color = INK
    return plot


def _add_markers(plot: figure, source: ColumnDataSource) -> Model:
    """Draw a marker source whose rows carry their own style."""
    return plot.scatter(
        "x",
        "y",
        marker="marker",
        size="size",
        angle="angle",
        fill_color="fill_color",
        fill_alpha="fill_alpha",
        line_color="line_color",
        line_width="line_width",
        source=source,
        level="annotation",
    )


def _add_labels(plot: figure, source: ColumnDataSource) -> None:
    """Draw a text source whose rows carry their own style."""
    plot.text(
        "x",
        "y",
        text="text",
        x_offset="x_offset",
        y_offset="y_offset",
        text_color="text_color",
        text_font_size="text_font_size",
        text_align="text_align",
        text_baseline="text_baseline",
        source=source,
        level="annotation",
    )


def _add_guides(plot: figure, source: ColumnDataSource) -> Model:
    """Draw vertical event guides shared by the linked signed-x panels."""
    return plot.vspan(
        x="x",
        line_color="line_color",
        line_dash="line_dash",
        line_alpha="line_alpha",
        line_width="line_width",
        source=source,
    )


def build_guided_plot(
    sources: dict[str, dict], properties: dict[str, dict]
) -> LayoutDOM:
    """Build the linked guided view filled with computed plot data.

    Parameters
    ----------
    sources, properties : dict
        Output of ``calculate_guided_plot_data``.

    Returns
    -------
    bokeh.models.LayoutDOM
        Projected velocity and gradient on the left; shell, vector
        construction, and radial law on the right.
    """
    source = {
        name: ColumnDataSource(data=data, name=name)
        for name, data in sources.items()
    }
    ranges = {
        name: Range1d(name=name)
        for name in (
            "projected-x",
            "projected-y",
            "gradient-y",
            "geometry-x",
            "geometry-y",
            "radial-x",
            "radial-y",
        )
    }

    projected = _figure(
        ranges["projected-x"],
        ranges["projected-y"],
        PROJECTED_WIDTH,
        PROJECTED_HEIGHT,
        title="Which forward crossings meet the line target?",
        min_border_left=70,
    )
    projected.xaxis.visible = False
    projected.yaxis.axis_label = "projected velocity  v∥(x)"
    projected.title.text_font_size = "12pt"
    projected.vstrip(
        x0="x0",
        x1="x1",
        fill_color="#e5f2ef",
        fill_alpha=0.75,
        line_color=None,
        source=source["active-shading"],
    )
    projected.hspan(y=0, line_color=MUTED, line_width=0.8)
    guides = _add_guides(projected, source["guides"])
    _add_guides(projected, source["guide-marks"])
    projected.hspan(
        y="y",
        line_color=RED,
        line_width=1.5,
        line_dash=DASHED,
        source=source["line-target"],
    )
    reference_line = projected.line(
        "x",
        "y",
        line_color=MUTED,
        line_width=1.6,
        line_dash=DASHED,
        line_alpha=0.75,
        source=source["projected-reference"],
    )
    projected.line(
        "x",
        "y",
        line_color=INK,
        line_width=1.5,
        line_alpha=0.65,
        source=source["projected-curve"],
    )
    active_line = projected.line(
        "x",
        "y",
        line_color=INK,
        line_width=3,
        line_cap="round",
        source=source["projected-active"],
    )
    markers = _add_markers(projected, source["projected-markers"])
    _add_markers(projected, source["projected-rings"])
    _add_labels(projected, source["projected-labels"])
    projected.add_layout(
        Label(
            x=0.98 * PROJECTED_WIDTH,
            y=0.97 * PROJECTED_HEIGHT,
            x_units="screen",
            y_units="screen",
            text_align="right",
            text_baseline="top",
            text_font_size="9pt",
            text_color=RED,
            name="line-target-label",
        )
    )
    projected.add_layout(
        Legend(
            items=[
                LegendItem(label={"field": "label"}, renderers=[renderer])
                for renderer in (active_line, reference_line, guides, markers)
            ],
            orientation="horizontal",
            ncols=4,
            label_text_font_size="8pt",
            label_text_color=INK,
            border_line_color=None,
            background_fill_alpha=0.0,
            glyph_width=16,
            spacing=8,
            margin=0,
        ),
        "above",
    )

    gradient = _figure(
        ranges["projected-x"],
        ranges["gradient-y"],
        PROJECTED_WIDTH,
        GRADIENT_HEIGHT,
        title="G(x) = dv∥/dx   ·   G>0: redshifting as packet travels; "
        "G<0: blueshifting as packet travels",
        min_border_left=70,
    )
    gradient.title.text_font_size = "8.5pt"
    gradient.yaxis.axis_label = "G(x)"
    gradient.hspan(y=0, line_color=MUTED, line_width=0.8)
    _add_guides(gradient, source["guides"])
    _add_guides(gradient, source["guide-marks"])
    gradient.line(
        "x",
        "y",
        line_color=MUTED,
        line_width=1.3,
        line_dash=DASHED,
        line_alpha=0.7,
        source=source["gradient-reference"],
    )
    gradient.line(
        "x",
        "y",
        line_color=GREEN,
        line_width=1.5,
        source=source["gradient-curve"],
    )
    gradient.line(
        "x",
        "y",
        line_color=GREEN,
        line_width=2.5,
        source=source["gradient-active"],
    )
    _add_markers(gradient, source["gradient-markers"])

    geometry = _figure(
        ranges["geometry-x"],
        ranges["geometry-y"],
        SIDE_WIDTH,
        SIDE_WIDTH,
        axes=False,
        title="Shell and ray",
        min_border_left=SIDE_BORDER,
    )
    geometry.annulus(
        x=0,
        y=0,
        inner_radius="inner_radius",
        outer_radius="outer_radius",
        fill_color="fill_color",
        line_color=None,
        source=source["shell-rings"],
    )
    _add_guides(geometry, source["guides"])
    _add_guides(geometry, source["guide-marks"])
    geometry.circle(
        x=0,
        y=0,
        radius="radius",
        fill_color="white",
        fill_alpha="fill_alpha",
        line_color=INK,
        line_width=1.2,
        source=source["shell-edges"],
    )
    geometry.circle(
        x=0,
        y=0,
        radius="radius",
        fill_alpha=0.0,
        line_color=MUTED,
        line_width=1.1,
        line_dash=SHORT_DASHED,
        source=source["stagnation-ring"],
    )
    geometry.hspan(
        y="y", line_color=MUTED, line_width=1, source=source["geometry-ray"]
    )
    geometry.segment(
        x0="x0",
        y0="y",
        x1="x1",
        y1="y",
        line_color=INK,
        line_width=2.7,
        source=source["geometry-active-ray"],
    )
    _add_markers(geometry, source["geometry-markers"])
    _add_labels(geometry, source["geometry-labels"])

    vector = _figure(
        ranges["geometry-x"],
        ranges["geometry-y"],
        SIDE_WIDTH,
        SIDE_WIDTH,
        axes=False,
        title="Velocity projection",
        min_border_left=SIDE_BORDER,
    )
    vector.segment(
        x0="x0",
        y0="y0",
        x1="x1",
        y1="y1",
        line_color="line_color",
        line_width="line_width",
        line_dash="line_dash",
        line_alpha="line_alpha",
        source=source["vector-segments"],
    )
    vector.add_layout(
        Arrow(
            x_start="x_start",
            y_start="y_start",
            x_end="x_end",
            y_end="y_end",
            line_color="line_color",
            line_width="line_width",
            end=NormalHead(
                size="head_size",
                fill_color="head_color",
                line_color="head_line_color",
                line_width=0.7,
            ),
            source=source["vector-arrows"],
            level="glyph",
        )
    )
    _add_markers(vector, source["vector-markers"])
    _add_labels(vector, source["vector-labels"])
    vector.add_layout(
        Label(
            x=0.5 * SIDE_WIDTH,
            y=0.62 * SIDE_WIDTH,
            x_units="screen",
            y_units="screen",
            text="No verified forward resonance on this flight",
            text_align="center",
            text_font_size="8pt",
            text_color=MUTED,
            name="vector-empty-note",
        )
    )
    vector.add_layout(
        Label(
            x=0.02 * SIDE_WIDTH,
            y=0.03 * SIDE_WIDTH,
            x_units="screen",
            y_units="screen",
            text_font_size="7pt",
            text_color=INK,
            name="vector-note",
        )
    )
    radial_key = vector.line(
        [np.nan], [np.nan], line_width=2, name="vector-legend-radial"
    )
    projection_key = vector.line(
        [np.nan], [np.nan], line_color=RED, line_width=2
    )
    vector.add_layout(
        Legend(
            items=[
                LegendItem(
                    label="radial vᵣ (shell color)", renderers=[radial_key]
                ),
                LegendItem(label="along ray v∥", renderers=[projection_key]),
            ],
            location="bottom_right",
            label_text_font_size="7pt",
            border_line_color=None,
            background_fill_alpha=0.0,
            glyph_width=14,
            label_height=10,
            spacing=0,
            padding=2,
        )
    )

    radial = _figure(
        ranges["radial-x"],
        ranges["radial-y"],
        SIDE_WIDTH,
        RADIAL_HEIGHT,
        title="v(r) = m r + q",
        x_axis_label="radius",
        y_axis_label="radial v",
        min_border_left=SIDE_BORDER,
    )
    radial.hspan(y=0, line_color=MUTED, line_width=0.8)
    radial_reference = radial.line(
        "x",
        "y",
        line_color=MUTED,
        line_width=1.5,
        line_dash=DASHED,
        line_alpha=0.75,
        source=source["radial-reference"],
    )
    radial_current = radial.line(
        "x",
        "y",
        line_color=INK,
        line_width=2.2,
        source=source["radial-current"],
    )
    radial.scatter("x", "y", size=6.5, color=INK, source=source["radial-edges"])
    radial.add_layout(
        Label(
            x=0.02 * SIDE_WIDTH,
            y=0.98 * RADIAL_HEIGHT,
            x_units="screen",
            y_units="screen",
            text_baseline="top",
            text_font_size="7.5pt",
            text_color=INK,
            name="radial-edges-label",
        )
    )
    radial.add_layout(
        Legend(
            items=[
                LegendItem(label={"field": "label"}, renderers=[renderer])
                for renderer in (radial_current, radial_reference)
            ],
            location="bottom_right",
            label_text_font_size="6.5pt",
            border_line_color=None,
            background_fill_alpha=0.0,
            glyph_width=14,
            label_height=10,
            spacing=0,
            padding=2,
        )
    )

    plot = row(column(projected, gradient), column(geometry, vector, radial))
    update_guided_plot(plot, {}, properties)
    return plot


def update_guided_plot(
    plot: LayoutDOM, sources: dict[str, dict], properties: dict[str, dict]
) -> None:
    """Refill a built guided view in place with new plot data.

    Parameters
    ----------
    plot : bokeh.models.LayoutDOM
        View returned by ``build_guided_plot``.
    sources, properties : dict
        Output of ``calculate_guided_plot_data``.
    """
    models = {model.name: model for model in plot.select({"type": Model})}
    for name, data in sources.items():
        models[name].data = data
    for name, values in properties.items():
        target = models[name]
        # The legend key line keeps its colour on the glyph, not the renderer
        (getattr(target, "glyph", target)).update(**values)


def render_guided_case(
    primary: CaseAnalysis,
    reference: CaseAnalysis | None = None,
    selected_slot_index: int | None = None,
    *,
    show_reference: bool = True,
) -> pn.pane.Bokeh:
    """Render the guided view for one flight as a notebook-ready pane.

    Parameters
    ----------
    primary : CaseAnalysis
        Current shell flight.
    reference : CaseAnalysis or None, optional
        Paired comparison drawn as dashed curves.
    selected_slot_index : int or None, optional
        Production slot to ring in all applicable views.
    show_reference : bool, optional
        Draw the paired comparison curves when a reference is supplied.

    Returns
    -------
    panel.pane.Bokeh
        Static linked view.
    """
    return pn.pane.Bokeh(
        build_guided_plot(
            *calculate_guided_plot_data(
                primary,
                reference,
                selected_slot_index,
                show_reference=show_reference,
            )
        )
    )
