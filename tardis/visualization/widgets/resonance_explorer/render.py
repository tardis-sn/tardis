"""Headless, linked views of one analyzed single-shell packet flight.

The renderer consumes an analysis snapshot. It samples the displayed equations
but does not classify roots, turns, or event ownership.
"""

from __future__ import annotations

from textwrap import fill
from typing import TYPE_CHECKING

import matplotlib as mpl
import numpy as np
import numpy.typing as npt
from matplotlib import patheffects
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.transforms import blended_transform_factory

if TYPE_CHECKING:
    from .analysis import CaseAnalysis

from .fingerprint import input_fingerprint

INK = "#30343b"
MUTED = "#899198"
RED = "#b44943"
GREEN = "#3b7667"
AMBER = "#a66b2e"
UNVERIFIED = "#a44d77"
STOP = "#74487c"
PALETTE = LinearSegmentedColormap.from_list(
    "shell_radial_velocity",
    mpl.colormaps["RdBu_r"](np.linspace(0.20, 0.80, 256)),
)


def _sample_flight(
    coordinates: npt.NDArray[np.float64], analysis: CaseAnalysis
) -> tuple[
    npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]
]:
    """Sample the equations shown in the linked curve panels."""
    impact = analysis.inputs.b
    radius = np.hypot(coordinates, impact)
    projected = np.full_like(coordinates, np.nan)
    gradient = np.full_like(coordinates, np.nan)
    np.divide(
        analysis.q * coordinates,
        radius,
        out=projected,
        where=radius > 0,
    )
    projected += analysis.m * coordinates
    np.divide(
        analysis.q * impact**2,
        radius**3,
        out=gradient,
        where=radius > 0,
    )
    gradient += analysis.m
    in_shell = (radius >= analysis.inputs.r_in) & (
        radius <= analysis.inputs.r_out
    )
    projected[~in_shell] = np.nan
    gradient[~in_shell] = np.nan
    return radius, projected, gradient


def _candidate_style(candidate: object) -> tuple[str, str, str]:
    """Choose marker shape, edge, and label from a candidate disposition."""
    if candidate.verified_event_eligible:
        return "o", RED, "verified eligible crossing"
    if candidate.accepted and not candidate.unsquared_verified:
        return "X", UNVERIFIED, "production accepted, unsquared unverified"
    if candidate.physical:
        return "x", AMBER, "physical root, event excluded"
    if candidate.x is not None and np.isfinite(candidate.x):
        return "^", MUTED, "squared-sign artifact"
    return "", MUTED, "unavailable slot"


def _draw_shell_field(
    ax: mpl.axes.Axes, analysis: CaseAnalysis, norm: TwoSlopeNorm
) -> None:
    """Draw the signed radial velocity field with concentric vector patches."""
    inputs = analysis.inputs
    for radius in np.linspace(inputs.r_out, inputs.r_in, 96):
        ax.add_patch(
            Circle(
                (0, 0),
                radius,
                facecolor=PALETTE(norm(analysis.m * radius + analysis.q)),
                edgecolor="none",
                zorder=0,
            )
        )
    ax.add_patch(
        Circle(
            (0, 0),
            inputs.r_in,
            facecolor="white",
            edgecolor=INK,
            lw=1.2,
            zorder=1,
        )
    )
    ax.add_patch(
        Circle(
            (0, 0),
            inputs.r_out,
            facecolor="none",
            edgecolor=INK,
            lw=1.2,
            zorder=1,
        )
    )
    stagnation = analysis.stagnation_radius
    if stagnation is not None and inputs.r_in <= stagnation <= inputs.r_out:
        ax.add_patch(
            Circle(
                (0, 0),
                stagnation,
                facecolor="none",
                edgecolor=MUTED,
                lw=1.1,
                ls=(0, (3, 2)),
                zorder=2,
            )
        )
        ax.text(
            0, -stagnation, "  radial v=0", va="center", fontsize=8, color=INK
        )


def _draw_event_guides(
    axes: tuple[mpl.axes.Axes, mpl.axes.Axes, mpl.axes.Axes],
    analysis: CaseAnalysis,
) -> None:
    """Repeat event coordinates in the shared signed-x axes."""
    for ax in axes:
        for coordinate, color, style in (
            (analysis.inputs.x_0, INK, (0, (3, 2))),
            (analysis.x_boundary, MUTED, (0, (3, 2))),
            (
                analysis.x_end,
                STOP if analysis.x_end != analysis.x_boundary else MUTED,
                (0, (2, 2)),
            ),
        ):
            ax.axvline(
                coordinate, color=color, lw=0.85, ls=style, alpha=0.65, zorder=1
            )
        for turn in analysis.turn_candidates:
            if np.isfinite(turn.x):
                turn_color = (
                    GREEN
                    if turn.transport_active
                    or (turn.active and "inflection" in turn.kind)
                    else AMBER
                    if turn.active
                    else MUTED
                )
                ax.axvline(
                    turn.x,
                    color=turn_color,
                    lw=0.65,
                    ls=(0, (1, 3)),
                    alpha=0.35,
                    zorder=1,
                )
        for crossing in analysis.crossings:
            if np.isfinite(crossing.x):
                ax.axvline(
                    crossing.x,
                    color=RED,
                    lw=0.75,
                    ls=(0, (1, 3)),
                    alpha=0.3,
                    zorder=1,
                )


def _draw_turns(
    axes: tuple[mpl.axes.Axes, mpl.axes.Axes, mpl.axes.Axes],
    analysis: CaseAnalysis,
) -> None:
    """Mark active and excluded gradient zeros in each applicable panel."""
    geometry_ax, projected_ax, gradient_ax = axes
    for turn_idx, turn in enumerate(analysis.turn_candidates):
        if not np.isfinite(turn.x):
            continue
        inflection = "inflection" in turn.kind.lower()
        marker = "s" if inflection else "D"
        color = (
            GREEN
            if turn.transport_active or (turn.active and inflection)
            else AMBER
            if turn.active
            else MUTED
        )
        face = color if turn.transport_active and not inflection else "white"
        for ax, ordinate in (
            (geometry_ax, analysis.inputs.b),
            (projected_ax, analysis.m * turn.x + analysis.q * turn.x / turn.r),
            (gradient_ax, 0.0),
        ):
            artist = ax.scatter(
                [turn.x],
                [ordinate],
                marker=marker,
                s=50,
                facecolors=face,
                edgecolors=color,
                linewidths=1.4,
                zorder=7,
            )
            artist.set_gid(f"turn-{turn_idx}")
        if turn in analysis.analytic_contacts:
            contact_marker = "P" if inflection else "*"
            for ax, ordinate in (
                (geometry_ax, analysis.inputs.b),
                (projected_ax, analysis.inputs.w),
            ):
                artist = ax.scatter(
                    [turn.x],
                    [ordinate],
                    marker=contact_marker,
                    s=135,
                    facecolors="white",
                    edgecolors=GREEN,
                    linewidths=1.6,
                    zorder=9,
                )
                artist.set_gid(f"target-contact-{turn_idx}")


def _draw_candidates(
    axes: tuple[mpl.axes.Axes, mpl.axes.Axes, mpl.axes.Axes],
    analysis: CaseAnalysis,
    selected_slot_index: int | None,
) -> None:
    """Place raw slots and numbered distinct crossings without regrouping."""
    geometry_ax, projected_ax, gradient_ax = axes
    crossing_coordinates = {crossing.x for crossing in analysis.crossings}
    slot_transform = blended_transform_factory(
        projected_ax.transData, projected_ax.transAxes
    )
    for candidate in analysis.candidates:
        if candidate.x is None or not np.isfinite(candidate.x):
            continue
        raw_marker = projected_ax.scatter(
            [candidate.x],
            [0.035],
            transform=slot_transform,
            marker="v",
            s=32,
            facecolors="white",
            edgecolors=MUTED,
            linewidths=1.1,
            zorder=7,
        )
        raw_marker.set_gid(f"raw-slot-{candidate.slot_idx}")
        marker, color, _ = _candidate_style(candidate)
        if not marker:
            continue
        # Numbered crossings own the prominent marker; raw-slot artifacts and
        # excluded roots retain their own shapes.
        if candidate.x not in crossing_coordinates or not candidate.physical:
            for ax, ordinate in (
                (geometry_ax, analysis.inputs.b),
                (projected_ax, candidate.projected_velocity),
            ):
                if ordinate is None or not np.isfinite(ordinate):
                    continue
                if marker == "x":
                    artist = ax.scatter(
                        [candidate.x],
                        [ordinate],
                        marker=marker,
                        s=43,
                        color=color,
                        linewidths=1.2,
                        zorder=8,
                    )
                else:
                    artist = ax.scatter(
                        [candidate.x],
                        [ordinate],
                        marker=marker,
                        s=43,
                        facecolors="white",
                        edgecolors=color,
                        linewidths=1.2,
                        zorder=8,
                    )
                artist.set_gid(f"candidate-slot-{candidate.slot_idx}")
        if candidate.slot_idx == selected_slot_index:
            for ax, ordinate in (
                (geometry_ax, analysis.inputs.b),
                (projected_ax, candidate.projected_velocity),
                (gradient_ax, candidate.gradient),
            ):
                if ordinate is None or not np.isfinite(ordinate):
                    continue
                artist = ax.scatter(
                    [candidate.x],
                    [ordinate],
                    marker="o",
                    s=180,
                    facecolors="none",
                    edgecolors=STOP,
                    linewidths=1.7,
                    zorder=10,
                )
                artist.set_gid(f"selected-slot-{candidate.slot_idx}")
    label_positions: list[tuple[float, int]] = []
    for number, crossing in enumerate(analysis.crossings, start=1):
        if crossing.x is None or not np.isfinite(crossing.x):
            continue
        pixel_x = geometry_ax.transData.transform(
            (crossing.x, analysis.inputs.b)
        )[0]
        label_tier = 0
        while any(
            tier == label_tier and abs(pixel_x - previous_x) < 15
            for previous_x, tier in label_positions
        ):
            label_tier += 1
        label_positions.append((pixel_x, label_tier))
        face = RED if crossing.verified_event_eligible else "white"
        for ax, ordinate in (
            (geometry_ax, analysis.inputs.b),
            (projected_ax, analysis.inputs.w),
        ):
            artist = ax.scatter(
                [crossing.x],
                [ordinate],
                marker="o",
                s=66,
                facecolors=face,
                edgecolors=RED,
                linewidths=1.4,
                zorder=9,
            )
            artist.set_gid(f"crossing-{number}")
            ax.annotate(
                str(number),
                (crossing.x, ordinate),
                xytext=(0, 9 + 13 * label_tier),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                color=RED,
                zorder=11,
            )
    nearest = analysis.nearest_resonance
    if nearest is not None and nearest.x is not None:
        for ax, ordinate in (
            (geometry_ax, analysis.inputs.b),
            (projected_ax, analysis.inputs.w),
        ):
            marker = ax.scatter(
                [nearest.x],
                [ordinate],
                marker="o",
                s=145,
                facecolors="none",
                edgecolors=INK,
                linewidths=1.7,
                zorder=9.5,
            )
            marker.set_gid("nearest-verified-line")


def draw_projection_construction(
    ax: mpl.axes.Axes,
    analysis: CaseAnalysis,
    *,
    compact: bool = False,
    axis_limits: tuple[tuple[float, float], tuple[float, float]] | None = None,
    color_norm: TwoSlopeNorm | None = None,
) -> None:
    """Draw radial velocities and their signed projections on the packet ray.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis reserved for the vector construction.
    analysis : CaseAnalysis
        Current flight. Only its distinct physical forward crossings are
        drawn; algebraic slots do not independently create vectors.
    compact : bool, optional
        Simplify labels for the small guided-view panel.
    axis_limits : tuple of tuple of float or None, optional
        Fixed x and y limits. The guided view passes the shell plot limits so
        moving the launch point cannot shift the plotted origin.
    color_norm : matplotlib.colors.TwoSlopeNorm or None, optional
        Shell-field color normalization. If omitted, use the current shell's
        edge speeds to construct the same normalization.
    """
    inputs = analysis.inputs
    if color_norm is None:
        speed_scale = max(abs(inputs.v_in), abs(inputs.v_out), 1e-9)
        color_norm = TwoSlopeNorm(
            vmin=-speed_scale, vcenter=0.0, vmax=speed_scale
        )
    impact = inputs.b
    chord_start, chord_end = analysis.outer_chord
    chord_span = max(chord_end - chord_start, 0.5)
    vector_data: list[tuple[float, float, float, float]] = []
    for crossing in analysis.crossings:
        if crossing.x is None:
            continue
        radius = np.hypot(crossing.x, impact)
        if radius <= 0:
            continue
        radial_speed = analysis.m * radius + analysis.q
        if abs(radial_speed) <= 1e-13 * max(
            1.0, abs(analysis.m * radius), abs(analysis.q)
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
        (0.32 if compact else 0.2) * chord_span / largest_speed
        if largest_speed > 0
        else float("inf")
    )
    vertical_scale = (
        (1.28 if compact else 0.8) * max(impact, 0.6) / largest_vertical
        if largest_vertical > 0
        else float("inf")
    )
    vector_scale = min(horizontal_scale, vertical_scale)
    if not np.isfinite(vector_scale):
        vector_scale = 1.0
    if impact == 0:
        for x, projected_speed, _, _ in vector_data:
            if projected_speed < 0 and x > inputs.x_0:
                vector_scale = min(
                    vector_scale, 0.7 * (x - inputs.x_0) / -projected_speed
                )
            elif projected_speed > 0 and x < analysis.x_end:
                vector_scale = min(
                    vector_scale, 0.7 * (analysis.x_end - x) / projected_speed
                )
    if axis_limits is not None:
        x_limits, y_limits = axis_limits
        x_margin = 0.03 * (x_limits[1] - x_limits[0])
        y_margin = 0.03 * (y_limits[1] - y_limits[0])
        for x, projected_speed, vertical_speed, _ in vector_data:
            for position, speed, limits, margin in (
                (x, projected_speed, x_limits, x_margin),
                (impact, vertical_speed, y_limits, y_margin),
            ):
                if speed > 0:
                    vector_scale = min(
                        vector_scale, (limits[1] - margin - position) / speed
                    )
                elif speed < 0:
                    vector_scale = min(
                        vector_scale, (position - limits[0] - margin) / -speed
                    )

    (full_ray,) = ax.plot(
        [chord_start, chord_end],
        [impact, impact],
        color=MUTED,
        lw=1.0,
        zorder=1,
    )
    full_ray.set_gid("vector-ray")
    (active_ray,) = ax.plot(
        [inputs.x_0, analysis.x_end],
        [impact, impact],
        color="#727a81",
        lw=1.1,
        solid_capstyle="round",
        zorder=2,
    )
    active_ray.set_gid("vector-active-trajectory")
    ax.scatter([0], [0], color=MUTED, s=16, zorder=2)
    ax.annotate(
        "O", (0, 0), xytext=(5, -12), textcoords="offset points", fontsize=9
    )
    for label, coordinate in (("start", inputs.x_0), ("end", analysis.x_end)):
        marker = ax.scatter(
            [coordinate],
            [impact],
            marker="|",
            color=MUTED,
            linewidths=1.0,
            s=85 if compact else 100,
            zorder=3,
        )
        marker.set_gid(f"vector-{label}")
    ax.annotate(
        "start",
        (inputs.x_0, impact),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=8 if compact else 9,
        color=MUTED,
        zorder=3,
    )
    ax.annotate(
        "end",
        (analysis.x_end, impact),
        xytext=(0, -16 if compact else 8),
        textcoords="offset points",
        ha="center",
        fontsize=8 if compact else 9,
        color=MUTED,
        zorder=3,
    )
    if analysis.x_end > inputs.x_0:
        arrow_start = inputs.x_0 + 0.70 * (analysis.x_end - inputs.x_0)
        arrow_end = inputs.x_0 + 0.76 * (analysis.x_end - inputs.x_0)
        ax.add_patch(
            FancyArrowPatch(
                (arrow_start, impact),
                (arrow_end, impact),
                arrowstyle="-|>",
                mutation_scale=10 if compact else 13,
                lw=1.1,
                color="#727a81",
                zorder=4,
            )
        )

    x_positions = (
        [inputs.x_0, analysis.x_end, 0.0]
        if compact
        else [chord_start, chord_end, 0.0]
    )
    y_positions = [0.0, impact]
    close_numbers: list[tuple[float, int]] = []
    clustered = False
    for number, (x, projected_speed, vertical_speed, radial_speed) in enumerate(
        vector_data, start=1
    ):
        radial_color = PALETTE(color_norm(radial_speed))
        point = (x, impact)
        radial_tip = (
            x + vector_scale * projected_speed,
            impact + vector_scale * vertical_speed,
        )
        projection_tip = (x + vector_scale * projected_speed, impact)
        x_positions.extend((radial_tip[0], projection_tip[0]))
        y_positions.append(radial_tip[1])
        ax.plot(
            [0, x],
            [0, impact],
            color=MUTED,
            lw=0.9,
            ls=(0, (2, 3)),
            alpha=0.55,
            zorder=1,
        )
        if radial_speed == 0.0:
            radial_marker = ax.scatter(
                [x],
                [impact],
                marker="o",
                s=85 if compact else 115,
                facecolors=radial_color,
                edgecolors="black",
                linewidths=0.7,
                zorder=6,
            )
            radial_marker.set_gid(f"vector-radial-{number}")
        else:
            radial_arrow = FancyArrowPatch(
                point,
                radial_tip,
                arrowstyle="-|>",
                mutation_scale=7 if compact else 11,
                shrinkA=0,
                shrinkB=0,
                lw=3.0 if impact == 0 else (1.8 if compact else 2.0),
                color=radial_color,
                zorder=6,
            )
            radial_arrow.set_gid(f"vector-radial-{number}")
            radial_arrow.set_path_effects(
                [
                    patheffects.withStroke(
                        linewidth=radial_arrow.get_linewidth() + 0.7,
                        foreground="black",
                    )
                ]
            )
            ax.add_patch(radial_arrow)
        if np.isclose(projected_speed, 0.0, rtol=0, atol=1e-14):
            projection = ax.scatter(
                [x],
                [impact],
                marker="o",
                s=90,
                facecolors="white",
                edgecolors=RED,
                linewidths=1.8,
                zorder=7,
            )
        else:
            projection = FancyArrowPatch(
                point,
                projection_tip,
                arrowstyle="-|>",
                mutation_scale=7 if compact else 11,
                shrinkA=0,
                shrinkB=0,
                lw=1.7 if compact else 1.9,
                color=RED,
                zorder=7,
            )
            ax.add_patch(projection)
        projection.set_gid(f"vector-projection-{number}")
        (perpendicular,) = ax.plot(
            [radial_tip[0], projection_tip[0]],
            [radial_tip[1], projection_tip[1]],
            color=MUTED,
            lw=1.0,
            ls=(0, (3, 2)),
            zorder=5,
        )
        perpendicular.set_gid(f"vector-perpendicular-{number}")
        marker = ax.scatter(
            [x],
            [impact],
            color=RED,
            edgecolors="white",
            linewidths=0.8,
            s=44 if compact else 72,
            zorder=9,
        )
        marker.set_gid(f"vector-point-{number}")
        label_tier = 0
        while any(
            tier == label_tier and abs(x - previous_x) < 0.07 * chord_span
            for previous_x, tier in close_numbers
        ):
            label_tier += 1
            clustered = True
        close_numbers.append((x, label_tier))
        ax.annotate(
            str(number),
            point,
            xytext=(0, 11 + 17 * label_tier),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.5 if compact else 9.5,
            color=RED,
            zorder=10,
        )
        if radial_speed == 0.0:
            ax.annotate(
                "vᵣ=0",
                point,
                xytext=(10, -14),
                textcoords="offset points",
                fontsize=7.5 if compact else 8.5,
                color=INK,
            )
        if projected_speed == 0.0:
            ax.annotate(
                "v∥=0",
                point,
                xytext=(10, 8),
                textcoords="offset points",
                fontsize=7.5 if compact else 8.5,
                color=RED,
            )

    if not vector_data:
        ax.text(
            0.5,
            0.62,
            "No verified forward resonance on this flight",
            transform=ax.transAxes,
            ha="center",
            fontsize=10,
            color=MUTED,
        )
    if compact:
        if impact == 0 and vector_data:
            ax.text(
                0.02,
                0.03,
                "vᵣ = v∥; arrows overlap",
                transform=ax.transAxes,
                fontsize=7,
                color=INK,
                va="bottom",
            )
        elif clustered:
            ax.text(
                0.02,
                0.03,
                "Nearby resonance labels",
                transform=ax.transAxes,
                fontsize=7,
                color=INK,
                va="bottom",
            )
    else:
        display_notes = ["All velocity arrows use one display scale."]
        if clustered:
            display_notes.append(
                "Nearby resonances: vectors overlap at this scale; numbers are staggered."
            )
        if impact == 0 and vector_data:
            display_notes.append(
                "vᵣ = v∥ on a radial ray (blue and red overlap)"
            )
        ax.text(
            0.02,
            0.04,
            "\n".join(display_notes),
            transform=ax.transAxes,
            fontsize=9,
            color=INK,
            va="bottom",
        )
    x_padding = max(0.08 * (max(x_positions) - min(x_positions)), 0.3)
    y_padding = max(0.10 * (max(y_positions) - min(y_positions)), 0.16)
    x_min = min(x_positions) - x_padding
    x_max = max(x_positions) + x_padding
    y_min = min(y_positions) - y_padding
    y_max = max(y_positions) + y_padding
    axis_box = ax.get_position()
    box_ratio = (
        axis_box.width
        * ax.figure.get_figwidth()
        / (axis_box.height * ax.figure.get_figheight())
    )
    x_span = x_max - x_min
    y_span = y_max - y_min
    if x_span / y_span < box_ratio:
        margin = (box_ratio * y_span - x_span) / 2
        x_min -= margin
        x_max += margin
    else:
        margin = (x_span / box_ratio - y_span) / 2
        y_min -= margin
        y_max += margin
    if axis_limits is None:
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
    else:
        ax.set_xlim(*axis_limits[0])
        ax.set_ylim(*axis_limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_yticks([])
    legend_radial_color = PALETTE(
        color_norm(vector_data[0][3] if vector_data else 0.0)
    )
    radial_legend = Line2D(
        [],
        [],
        color=legend_radial_color,
        lw=2,
        label="radial vᵣ (shell color)",
    )
    radial_legend.set_path_effects(
        [patheffects.withStroke(linewidth=2.7, foreground="black")]
    )
    if compact:
        ax.set_xticks([])
        ax.set_title(
            "Velocity projection", loc="left", fontsize=9, y=0.96, pad=0
        )
        ax.legend(
            handles=[
                radial_legend,
                Line2D([], [], color=RED, lw=2, label="along ray v∥"),
            ],
            loc="lower right",
            fontsize=7,
            frameon=False,
            handlelength=1.4,
            labelspacing=0.2,
            borderpad=0.1,
        )
    else:
        ax.set_xlabel("signed x   ·   packet moves right →", fontsize=10)
        ax.set_title(
            f"Radial velocity and its projection at each resonance   ·   v∥ = vᵣ x/r = w = {inputs.w:.5g}",
            loc="left",
            fontsize=11,
        )
        ax.legend(
            handles=[
                radial_legend,
                Line2D([], [], color=RED, lw=2.4, label="ray projection v∥"),
            ],
            loc="lower right",
            fontsize=9,
            frameon=False,
        )


def render_guided_case(
    primary: CaseAnalysis,
    reference: CaseAnalysis | None = None,
    selected_slot_index: int | None = None,
    *,
    show_reference: bool = True,
) -> Figure:
    """Render a compact, projection-led lesson view.

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
    Figure
        Agg-backed figure with a dominant projected-velocity view.
    """
    inputs = primary.inputs
    cases = (primary,) if reference is None else (primary, reference)
    outer_edge = max(
        np.sqrt(max(case.inputs.r_out**2 - case.inputs.b**2, 0.0))
        for case in cases
    )
    half_width = max(1.08 * max(case.inputs.r_out for case in cases), 0.5)
    x_limits = (-half_width, half_width)
    coordinates = np.linspace(-outer_edge, outer_edge, 1201)
    sampled = [_sample_flight(coordinates, case) for case in cases]
    radial_values = [
        case.m * np.linspace(case.inputs.r_in, case.inputs.r_out, 220) + case.q
        for case in cases
    ]
    projected_extrema = [inputs.w, 0.0]
    gradient_extrema = [0.0]
    radial_extrema = [0.0]
    for (_, projected, gradient), radial in zip(
        sampled, radial_values, strict=True
    ):
        projected_extrema.extend(projected[np.isfinite(projected)])
        gradient_extrema.extend(gradient[np.isfinite(gradient)])
        radial_extrema.extend(radial)

    def padded_limits(
        values: list[float] | npt.NDArray[np.float64],
    ) -> tuple[float, float]:
        lower = float(np.min(values))
        upper = float(np.max(values))
        padding = max(0.08 * (upper - lower), 0.12)
        return lower - padding, upper + padding

    figure = Figure(figsize=(12.8, 8.4), dpi=130, facecolor="white")
    FigureCanvasAgg(figure)
    projected_ax = figure.add_axes((0.07, 0.32, 0.64, 0.62))
    gradient_ax = figure.add_axes((0.07, 0.09, 0.64, 0.16), sharex=projected_ax)
    geometry_ax = figure.add_axes((0.75, 0.64, 0.22, 0.30))
    vector_ax = figure.add_axes((0.75, 0.31, 0.22, 0.30))
    radial_ax = figure.add_axes((0.75, 0.09, 0.22, 0.19))
    for ax, name in (
        (projected_ax, "guided projected velocity"),
        (gradient_ax, "guided projected gradient"),
        (geometry_ax, "guided shell geometry"),
        (radial_ax, "guided radial law"),
        (vector_ax, "guided vector construction"),
    ):
        ax.set_gid(name)
        ax.tick_params(colors=INK, labelsize=8)
        for spine in ax.spines.values():
            spine.set_color(MUTED)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    projected_ax.spines["bottom"].set_visible(False)
    for ax in (geometry_ax, vector_ax):
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    active_shading = projected_ax.axvspan(
        inputs.x_0, primary.x_end, facecolor="#e5f2ef", alpha=0.75, zorder=0
    )
    active_shading.set_gid("active-flight-shading")
    target = projected_ax.axhline(
        inputs.w, color=RED, lw=1.5, ls="--", zorder=2
    )
    target.set_gid("line-target")
    projected_ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    gradient_ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    projected_ax.set_xlim(*x_limits)
    projected_ax.set_ylim(*padded_limits(projected_extrema))
    gradient_ax.set_ylim(*padded_limits(gradient_extrema))
    projected_ax.set_ylabel("projected velocity  v∥(x)", fontsize=9)
    projected_ax.set_title(
        "Which forward crossings meet the line target?", loc="left", fontsize=12
    )
    projected_ax.annotate(
        f"line target w = {inputs.w:.6g}",
        (0.98, 0.97),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=9,
        color=RED,
    )
    gradient_ax.set_ylabel("G(x)", fontsize=9)
    gradient_ax.set_title(
        "G(x) = dv∥/dx   ·   G>0: redshifting as packet travels; G<0: blueshifting as packet travels",
        loc="left",
        fontsize=8.5,
    )
    gradient_ax.set_xlabel("")
    projected_ax.tick_params(
        axis="x", which="both", bottom=False, labelbottom=False
    )

    current_label = "homology" if primary.q == 0 else "nonhomology"
    if reference is not None and show_reference:
        reference_label = "homology" if reference.q == 0 else "nonhomology"
        projected_ax.plot(
            coordinates,
            sampled[1][1],
            color=MUTED,
            lw=1.6,
            ls=(0, (4, 3)),
            alpha=0.75,
            label=f"{reference_label} reference",
            zorder=2,
        )
        gradient_ax.plot(
            coordinates,
            sampled[1][2],
            color=MUTED,
            lw=1.3,
            ls=(0, (4, 3)),
            alpha=0.7,
            zorder=2,
        )
    projected_ax.plot(
        coordinates, sampled[0][1], color=INK, lw=1.5, alpha=0.65, zorder=3
    )
    gradient_ax.plot(coordinates, sampled[0][2], color=GREEN, lw=1.5, zorder=3)
    active_coordinates = np.linspace(inputs.x_0, primary.x_end, 500)
    _, active_projection, active_gradient = _sample_flight(
        active_coordinates, primary
    )
    (active_projected,) = projected_ax.plot(
        active_coordinates,
        active_projection,
        color=INK,
        lw=3,
        solid_capstyle="round",
        zorder=5,
        label=f"current {current_label}",
    )
    active_projected.set_gid("active-projection")
    (active_gradient_artist,) = gradient_ax.plot(
        active_coordinates, active_gradient, color=GREEN, lw=2.5, zorder=4
    )
    active_gradient_artist.set_gid("active-gradient")
    if reference is not None and show_reference:
        projected_ax.legend(loc="upper left", fontsize=8, frameon=False)

    speed_scale = max(abs(value) for value in radial_extrema) or 1e-9
    color_norm = TwoSlopeNorm(vmin=-speed_scale, vcenter=0, vmax=speed_scale)
    _draw_shell_field(geometry_ax, primary, color_norm)
    geometry_ax.set_xlim(*x_limits)
    geometry_ax.set_ylim(-half_width, half_width)
    geometry_ax.set_aspect("equal", adjustable="box")
    geometry_ax.set_title("Shell and ray", loc="left", fontsize=10)
    geometry_ax.set_xlabel("")
    geometry_ax.tick_params(
        axis="x", which="both", bottom=False, labelbottom=False
    )
    geometry_ax.set_yticks([])
    geometry_ax.axhline(inputs.b, color=MUTED, lw=1, zorder=2)
    (active_ray,) = geometry_ax.plot(
        [inputs.x_0, primary.x_end],
        [inputs.b] * 2,
        color=INK,
        lw=2.7,
        zorder=4,
    )
    active_ray.set_gid("active-ray")
    geometry_ax.scatter([inputs.x_0], [inputs.b], color=INK, s=34, zorder=7)
    geometry_ax.scatter(
        [primary.x_end],
        [inputs.b],
        facecolors="white",
        edgecolors=INK,
        marker="s",
        s=37,
        zorder=7,
    )
    geometry_ax.annotate(
        "start",
        (inputs.x_0, inputs.b),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=7,
    )
    geometry_ax.annotate(
        "end",
        (primary.x_end, inputs.b),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=7,
    )

    radii = np.linspace(inputs.r_in, inputs.r_out, 220)
    radial_ax.axhline(0, color=MUTED, lw=0.8)
    (radial_line,) = radial_ax.plot(
        radii,
        radial_values[0],
        color=INK,
        lw=2.2,
        zorder=4,
        label=f"current {current_label}",
    )
    radial_line.set_gid("current-radial-law")
    if reference is not None and show_reference:
        radial_ax.plot(
            radii,
            radial_values[1],
            color=MUTED,
            lw=1.5,
            ls=(0, (4, 3)),
            alpha=0.75,
            label=f"{reference_label} reference",
        )
        radial_ax.legend(loc="lower right", fontsize=6.5, frameon=False)
    radial_ax.scatter(
        [inputs.r_in, inputs.r_out],
        [inputs.v_in, inputs.v_out],
        color=INK,
        s=24,
        zorder=5,
    )
    radial_ax.set_xlim(inputs.r_in - 0.08, inputs.r_out + 0.08)
    radial_ax.set_ylim(*padded_limits(radial_extrema))
    radial_ax.set_title("v(r) = m r + q", loc="left", fontsize=10)
    radial_ax.set_xlabel("radius", fontsize=8)
    radial_ax.set_ylabel("radial v", fontsize=8)
    radial_ax.text(
        0.02,
        0.98,
        f"edges: {inputs.v_in:.5g}, {inputs.v_out:.5g}",
        transform=radial_ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color=INK,
    )

    linked_axes = (geometry_ax, projected_ax, gradient_ax)
    _draw_event_guides(linked_axes, primary)
    _draw_turns(linked_axes, primary)
    _draw_candidates(linked_axes, primary, selected_slot_index)
    draw_projection_construction(
        vector_ax,
        primary,
        compact=True,
        axis_limits=(geometry_ax.get_xlim(), geometry_ax.get_ylim()),
        color_norm=color_norm,
    )
    return figure


def render_case(
    analysis: CaseAnalysis, selected_slot_index: int | None = None
) -> Figure:
    """Render a complete analysis snapshot as a headless linked figure.

    Parameters
    ----------
    analysis : CaseAnalysis
        One validated and immutable shell-flight analysis.
    selected_slot_index : int or None, optional
        Raw production slot to ring in each applicable view.

    Returns
    -------
    Figure
        Agg-backed figure ready for Panel display or PDF/PNG export.
    """
    inputs = analysis.inputs
    edge = float(np.sqrt(max(inputs.r_out**2 - inputs.b**2, 0.0)))
    half_width = max(1.1 * inputs.r_out, 0.5)
    x_limits = (-half_width, half_width)
    speed_scale = max(abs(inputs.v_in), abs(inputs.v_out), 1e-9)
    color_norm = TwoSlopeNorm(vmin=-speed_scale, vcenter=0.0, vmax=speed_scale)
    figure = Figure(figsize=(13.3, 20.0), dpi=130, facecolor="white")
    FigureCanvasAgg(figure)
    geometry_ax = figure.add_axes((0.08, 0.48, 0.62, 0.48))
    radial_ax = figure.add_axes((0.76, 0.75, 0.20, 0.19))
    projected_ax = figure.add_axes((0.08, 0.31, 0.62, 0.13), sharex=geometry_ax)
    gradient_ax = figure.add_axes((0.08, 0.20, 0.62, 0.08), sharex=geometry_ax)
    vector_ax = figure.add_axes((0.08, 0.04, 0.88, 0.12))
    for ax, name in (
        (geometry_ax, "shell geometry"),
        (radial_ax, "radial law"),
        (projected_ax, "projected velocity"),
        (gradient_ax, "projected gradient"),
        (vector_ax, "vector construction"),
    ):
        ax.set_gid(name)
        ax.tick_params(colors=INK, labelsize=8)
        for spine in ax.spines.values():
            spine.set_color(MUTED)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    projected_ax.spines["bottom"].set_visible(False)
    geometry_ax.spines["left"].set_visible(False)
    geometry_ax.spines["bottom"].set_visible(False)
    vector_ax.spines["left"].set_visible(False)

    _draw_shell_field(geometry_ax, analysis, color_norm)
    geometry_ax.set_xlim(*x_limits)
    geometry_ax.set_ylim(-inputs.r_out * 1.12, inputs.r_out * 1.12)
    geometry_ax.set_aspect("equal", adjustable="box")
    geometry_ax.set_yticks([])
    geometry_ax.tick_params(
        axis="x", which="both", bottom=False, labelbottom=False
    )
    geometry_ax.set_title(
        "One shell and the packet ray", loc="left", fontsize=11
    )
    geometry_ax.axhline(inputs.b, color=MUTED, lw=1.5, zorder=3)
    (active_ray,) = geometry_ax.plot(
        [inputs.x_0, analysis.x_end],
        [inputs.b] * 2,
        color=INK,
        lw=3.1,
        solid_capstyle="round",
        zorder=4,
    )
    active_ray.set_gid("active-ray")
    if analysis.x_end > inputs.x_0:
        arrow_x = inputs.x_0 + 0.58 * (analysis.x_end - inputs.x_0)
        arrow_length = min(
            0.10 * (analysis.x_end - inputs.x_0), 0.12 * inputs.r_out
        )
        geometry_ax.add_patch(
            FancyArrowPatch(
                (arrow_x - arrow_length, inputs.b),
                (arrow_x + arrow_length, inputs.b),
                arrowstyle="-|>",
                mutation_scale=15,
                lw=1.6,
                color=INK,
                zorder=6,
            )
        )
    geometry_ax.scatter(
        [0], [inputs.b], marker="|", s=90, color=MUTED, zorder=5
    )
    geometry_ax.annotate(
        "closest approach",
        (0, inputs.b),
        xytext=(4, -18),
        textcoords="offset points",
        fontsize=7.5,
        color=INK,
    )
    geometry_ax.scatter(
        [inputs.x_0], [inputs.b], marker="o", s=68, color=INK, zorder=8
    )
    geometry_ax.annotate(
        "start",
        (inputs.x_0, inputs.b),
        xytext=(0, 10),
        textcoords="offset points",
        ha="center",
        fontsize=8,
    )
    geometry_ax.scatter(
        [analysis.x_boundary],
        [inputs.b],
        marker="s",
        s=64,
        facecolors="white",
        edgecolors=INK,
        zorder=7,
    )
    geometry_ax.annotate(
        f"{analysis.boundary_kind} boundary",
        (analysis.x_boundary, inputs.b),
        xytext=(0, -23),
        textcoords="offset points",
        ha="center",
        fontsize=7.7,
    )
    if analysis.x_end != analysis.x_boundary:
        geometry_ax.scatter(
            [analysis.x_end], [inputs.b], marker="8", s=82, color=STOP, zorder=8
        )
        geometry_ax.annotate(
            analysis.stop_reason,
            (analysis.x_end, inputs.b),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=STOP,
        )

    radii = np.linspace(inputs.r_in, inputs.r_out, 220)
    radial_speed = analysis.m * radii + analysis.q
    radial_ax.axhline(0, color=MUTED, lw=0.9)
    radial_ax.plot(radii, radial_speed, color=INK, lw=2.2)
    radial_ax.scatter(
        [inputs.r_in, inputs.r_out],
        [inputs.v_in, inputs.v_out],
        color=INK,
        s=21,
        zorder=4,
    )
    radial_ax.set_xlabel("radius", fontsize=9)
    radial_ax.set_ylabel("radial v(r)", fontsize=9)
    radial_ax.set_title("v(r) = m r + q", fontsize=10)
    if (
        analysis.stagnation_radius is not None
        and inputs.r_in <= analysis.stagnation_radius <= inputs.r_out
    ):
        radial_ax.scatter(
            [analysis.stagnation_radius],
            [0],
            marker="o",
            s=48,
            facecolors="white",
            edgecolors=INK,
            zorder=5,
        )
        radial_ax.annotate(
            "v=0",
            (analysis.stagnation_radius, 0),
            xytext=(5, 6),
            textcoords="offset points",
            fontsize=8,
        )
    sm = mpl.cm.ScalarMappable(norm=color_norm, cmap=PALETTE)
    colorbar = figure.colorbar(sm, ax=radial_ax, fraction=0.06, pad=0.08)
    colorbar.set_label("radial velocity", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

    coordinates = np.linspace(-edge, edge, 1201)
    _, projected_values, gradient_values = _sample_flight(coordinates, analysis)
    projected_ax.axhline(
        inputs.w,
        color=RED,
        lw=1.4,
        ls="--",
        label=f"line target w = {inputs.w:.6g}",
    )
    projected_ax.axhline(0, color=MUTED, lw=0.8)
    projected_ax.plot(coordinates, projected_values, color=INK, lw=1.9)
    active_coordinates = np.linspace(inputs.x_0, analysis.x_end, 500)
    _, active_projection, active_gradient = _sample_flight(
        active_coordinates, analysis
    )
    (active_projected,) = projected_ax.plot(
        active_coordinates, active_projection, color=INK, lw=3.0
    )
    active_projected.set_gid("active-projection")
    projected_ax.set_ylabel("projected v∥(x)", fontsize=9)
    projected_ax.set_title(
        "Projected velocity and line target", loc="left", fontsize=10
    )
    projected_ax.legend(loc="best", fontsize=8, frameon=False)
    gradient_ax.axhline(0, color=MUTED, lw=0.9)
    gradient_ax.plot(coordinates, gradient_values, color=GREEN, lw=1.5)
    (active_gradient_artist,) = gradient_ax.plot(
        active_coordinates, active_gradient, color=GREEN, lw=2.8
    )
    active_gradient_artist.set_gid("active-gradient")
    gradient_ax.set_ylabel("G(x) = dv∥/dx", fontsize=9)
    gradient_ax.set_xlabel(
        "signed x   (forward packet distance s = x − x₀)", fontsize=9
    )
    gradient_ax.set_title("Projected gradient", loc="left", fontsize=10)
    gradient_ax.set_xlim(*x_limits)

    linked_axes = (geometry_ax, projected_ax, gradient_ax)
    _draw_event_guides(linked_axes, analysis)
    _draw_turns(linked_axes, analysis)
    _draw_candidates(linked_axes, analysis, selected_slot_index)
    draw_projection_construction(vector_ax, analysis, color_norm=color_norm)
    projected_ax.tick_params(
        axis="x", which="both", bottom=False, labelbottom=False
    )
    legend_handles = [
        Line2D(
            [],
            [],
            color=INK,
            marker="o",
            ls="",
            markersize=7,
            label="packet start",
        ),
        Line2D(
            [],
            [],
            color=GREEN,
            marker="D",
            markerfacecolor=GREEN,
            ls="",
            markersize=6,
            label="transport interval reversal",
        ),
        Line2D(
            [],
            [],
            color=AMBER,
            marker="D",
            markerfacecolor="white",
            ls="",
            markersize=6,
            label="geometric turn, transport excluded",
        ),
        Line2D(
            [],
            [],
            color=MUTED,
            marker="D",
            markerfacecolor="white",
            ls="",
            markersize=6,
            label="turn outside active flight",
        ),
        Line2D(
            [],
            [],
            color=GREEN,
            marker="s",
            markerfacecolor="white",
            ls="",
            markersize=6,
            label="stationary inflection",
        ),
        Line2D(
            [],
            [],
            color=GREEN,
            marker="*",
            markerfacecolor="white",
            ls="",
            markersize=9,
            label="target at extremum",
        ),
        Line2D(
            [],
            [],
            color=MUTED,
            marker="v",
            markerfacecolor="white",
            ls="",
            markersize=6,
            label="raw solver slot (strip)",
        ),
        Line2D(
            [],
            [],
            color=RED,
            marker="o",
            ls="",
            markersize=7,
            label="numbered verified crossing",
        ),
        Line2D(
            [],
            [],
            color=INK,
            marker="o",
            markerfacecolor="none",
            ls="",
            markersize=9,
            label="nearest verified line",
        ),
        Line2D(
            [],
            [],
            color=UNVERIFIED,
            marker="X",
            ls="",
            markersize=7,
            label="production accepted, unsquared unverified",
        ),
        Line2D(
            [],
            [],
            color=AMBER,
            marker="x",
            ls="",
            markersize=7,
            label="physical, event excluded",
        ),
        Line2D(
            [],
            [],
            color=MUTED,
            marker="^",
            markerfacecolor="white",
            ls="",
            markersize=7,
            label="squared-sign artifact",
        ),
        Line2D(
            [],
            [],
            color=STOP,
            marker="o",
            markerfacecolor="none",
            ls="",
            markersize=9,
            label="selected raw slot",
        ),
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.76, 0.67),
        frameon=False,
        fontsize=8.2,
    )
    nearest = analysis.nearest_resonance
    production_nearest = analysis.production_nearest_resonance
    production_line = (
        "production next: none"
        if production_nearest is None
        else f"production next: slot {production_nearest.slot_idx} "
        f"({'verified' if production_nearest.unsquared_verified else 'unverified'})"
    )
    nearest_line = (
        "nearest verified resonance: none"
        if nearest is None
        else f"nearest verified resonance: s={nearest.s:.5g}"
    )
    status = "\n".join(
        (
            fill(
                f"{analysis.ray_family} → {analysis.shape_family} → {analysis.result_family}",
                width=28,
            ),
            fill(
                f"{analysis.radial_law}; solver: {analysis.solver_route}",
                width=28,
            ),
            f"shell: [{inputs.r_in:.5g}, {inputs.r_out:.5g}]",
            f"edge speeds: {inputs.v_in:.5g}, {inputs.v_out:.5g}",
            f"b={inputs.b:.5g}; x₀={inputs.x_0:.5g}; w={inputs.w:.5g}",
            f"active x: ({inputs.x_0:.5g}, {analysis.x_end:.5g}]",
            f"s_end={analysis.active_end:.5g}",
            f"first boundary: {analysis.boundary_kind}",
            f"geometric turns: {len(analysis.turns)}; packet splits: "
            f"{len(analysis.transport_turns)}",
            f"verified crossings: {len(analysis.crossings)}",
            fill(production_line, width=28),
            fill(nearest_line, width=28),
        )
    )
    figure.text(
        0.76,
        0.48,
        status,
        va="top",
        fontsize=8.1,
        color=INK,
        linespacing=1.5,
        wrap=True,
    )
    figure.text(
        0.76,
        0.30,
        "Filled red: verified eligible line.\n"
        "Open red: crossing excluded\n"
        "by event ownership.\n"
        "A resonance is not a sampled\n"
        "absorption. Shell color shows\n"
        "radial velocity.",
        va="top",
        fontsize=7.7,
        color=INK,
        linespacing=1.45,
    )
    figure.text(
        0.08,
        0.025,
        f"input fingerprint: {input_fingerprint(inputs)}",
        fontsize=7.5,
        color=MUTED,
    )
    figure.suptitle(
        "Single-shell nonhomology resonance explorer",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=13,
        color=INK,
    )
    return figure
