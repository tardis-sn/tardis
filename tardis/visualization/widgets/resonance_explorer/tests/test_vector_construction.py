"""Scientific and display checks for the radial-vector construction view."""

import math

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch

from tardis.visualization.widgets.resonance_explorer import (
    ExplorerInputs,
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer import render as rendering


def linked_inputs(**changes: float) -> ExplorerInputs:
    """Return the three-crossing oblique ray in illustrative units."""
    values = dict(
        r_in=1.0,
        r_out=3.0,
        v_in=-1.34375,
        v_out=0.65625,
        b=1.2,
        x_0=-math.sqrt(7.56),
        w=0.3,
    )
    values.update(changes)
    return ExplorerInputs(**values)


def vector_axis(case: object) -> object:
    """Draw only the vector construction for direct geometric assertions."""
    figure = Figure(figsize=(10, 3), dpi=100)
    FigureCanvasAgg(figure)
    axis = figure.add_subplot()
    rendering.draw_projection_construction(axis, case)
    figure.canvas.draw()
    return axis


def artist(axis: object, gid: str) -> object:
    """Find a scientific vector primitive by its stable identifier."""
    return next(item for item in axis.get_children() if item.get_gid() == gid)


def arrow_endpoints(axis: object, gid: str) -> tuple[np.ndarray, np.ndarray]:
    """Read the mathematical endpoints, excluding the drawn arrowhead."""
    arrow = artist(axis, gid)
    assert isinstance(arrow, FancyArrowPatch)
    start, end = arrow._posA_posB
    return np.asarray(start, dtype=float), np.asarray(end, dtype=float)


def test_linked_vector_tips_obey_independent_radial_projection() -> None:
    """Claim: each radial vector projects onto the directed ray by v_r x/r.

    Regime: three distinct forward crossings, including inbound negative flow.
    Verification: radial and projected arrow endpoints obey the dot product
    of v_r(x,b)/r with the fixed packet direction (1,0). A single positive
    velocity-to-length scale applies at every crossing. Exact artist endpoints
    are compared at roundoff precision; no sampled curve values are used.
    """
    case = analyze_case(linked_inputs())
    assert len(case.crossings) == 3
    axis = vector_axis(case)
    scales = []
    for number, crossing in enumerate(case.crossings, start=1):
        x = crossing.x
        assert x is not None
        point = np.array([x, case.inputs.b])
        radius = math.hypot(x, case.inputs.b)
        radial_speed = case.m * radius + case.q
        radial_start, radial_end = arrow_endpoints(
            axis, f"vector-radial-{number}"
        )
        projected_start, projected_end = arrow_endpoints(
            axis, f"vector-projection-{number}"
        )
        np.testing.assert_allclose(radial_start, point, rtol=0, atol=1e-12)
        np.testing.assert_allclose(projected_start, point, rtol=0, atol=1e-12)
        radial_delta = radial_end - radial_start
        projected_delta = projected_end - projected_start
        scale = np.linalg.norm(radial_delta) / abs(radial_speed)
        scales.append(scale)
        np.testing.assert_allclose(
            radial_delta,
            scale * radial_speed * np.array([x, case.inputs.b]) / radius,
            rtol=1e-12,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            projected_delta,
            [scale * radial_speed * x / radius, 0],
            rtol=1e-12,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            projected_delta[0], radial_delta[0], atol=1e-12
        )
        np.testing.assert_allclose(
            artist(axis, f"vector-point-{number}").get_offsets()[0], point
        )
        artist(axis, f"vector-perpendicular-{number}")
    assert min(scales) > 0
    np.testing.assert_allclose(scales, scales[0], rtol=1e-12)


def test_collapse_radial_arrows_point_inward_while_projection_points_forward() -> (
    None
):
    """Claim: negative radial speed does not imply negative ray projection.

    Regime: whole-shell collapse and two inbound crossings. Verification:
    v_r<0 at both crossings; the inward radial arrows point right/down while
    their signed ray projections point right toward increasing x.
    """
    case = analyze_case(linked_inputs(v_in=-2.5, v_out=-0.5, w=0.6))
    assert len(case.crossings) == 2
    axis = vector_axis(case)
    for number in (1, 2):
        radial_start, radial_end = arrow_endpoints(
            axis, f"vector-radial-{number}"
        )
        projected_start, projected_end = arrow_endpoints(
            axis, f"vector-projection-{number}"
        )
        assert radial_end[0] > radial_start[0]
        assert radial_end[1] < radial_start[1]
        assert projected_end[0] > projected_start[0]


def test_negative_line_target_reverses_the_projection_arrow() -> None:
    """Claim: the projected arrow retains the sign of a negative line target.

    Regime: homologous expansion, m=1, q=0, w=-0.3; exact root x=w.
    Verification: dot(v_r(x,b)/r,(1,0))=x=-0.3 points against packet travel.
    """
    case = analyze_case(linked_inputs(v_in=1, v_out=3, w=-0.3))
    assert len(case.crossings) == 1
    np.testing.assert_allclose(case.crossings[0].x, -0.3, atol=1e-14)
    axis = vector_axis(case)
    start, end = arrow_endpoints(axis, "vector-projection-1")
    assert end[0] < start[0]
    np.testing.assert_allclose(end[1], start[1], atol=1e-12)


def test_zero_projection_is_marked_without_a_direction_arrow() -> None:
    """Claim: at closest approach, radial flow has zero ray projection.

    Regime: homologous expansion with w=0 and b>0. The exact root is x=0;
    a point marker expresses the zero vector without implying a direction.
    """
    case = analyze_case(linked_inputs(v_in=1, v_out=3, w=0.0))
    assert len(case.crossings) == 1
    np.testing.assert_allclose(case.crossings[0].x, 0.0, atol=1e-14)
    axis = vector_axis(case)
    radial_start, radial_end = arrow_endpoints(axis, "vector-radial-1")
    np.testing.assert_allclose(radial_end[0], radial_start[0], atol=1e-12)
    assert radial_end[1] > radial_start[1]
    projection = artist(axis, "vector-projection-1")
    assert not isinstance(projection, FancyArrowPatch)
    np.testing.assert_allclose(
        projection.get_offsets()[0], [0, case.inputs.b], atol=1e-12
    )


def test_stagnation_resonances_have_no_spurious_vector_direction() -> None:
    """Claim: zero radial speed means both velocity vectors vanish.

    Regime: the nonhomologous linked law at w=0, with exact stagnation radius
    r=-q/m=2.34375. The inbound and outbound roots at this radius have
    v_r=v_parallel=0. Verification: their arrows collapse to a point or are
    drawn as zero markers; neither may point in an invented direction.
    """
    case = analyze_case(linked_inputs(w=0.0))
    assert len(case.crossings) == 3
    axis = vector_axis(case)
    for number in (1, 3):
        point = np.array([case.crossings[number - 1].x, case.inputs.b])
        for prefix in ("vector-radial", "vector-projection"):
            vector = artist(axis, f"{prefix}-{number}")
            if isinstance(vector, FancyArrowPatch):
                start, end = arrow_endpoints(axis, f"{prefix}-{number}")
                np.testing.assert_allclose(start, point, atol=1e-12)
                np.testing.assert_allclose(end, point, atol=1e-12)
            else:
                np.testing.assert_allclose(vector.get_offsets()[0], point)


def test_radial_ray_has_coincident_radial_and_projected_vectors() -> None:
    """Claim: b=0 makes the radial velocity entirely parallel to packet travel.

    Regime: homologous expansion, inward radial ray from x=-3 to the inner
    boundary, exact resonance x=-2 and w=-2. Verification: both vector tips
    are at the same signed-ray location, independent of display scaling.
    """
    case = analyze_case(ExplorerInputs(1, 3, 1, 3, 0, -3, -2))
    assert len(case.crossings) == 1
    axis = vector_axis(case)
    radial_start, radial_end = arrow_endpoints(axis, "vector-radial-1")
    projected_start, projected_end = arrow_endpoints(
        axis, "vector-projection-1"
    )
    np.testing.assert_allclose(radial_start, [-2, 0], atol=1e-12)
    np.testing.assert_allclose(projected_start, radial_start, atol=1e-12)
    np.testing.assert_allclose(projected_end, radial_end, atol=1e-12)
    assert projected_end[0] > case.inputs.x_0
    radial = artist(axis, "vector-radial-1")
    projection = artist(axis, "vector-projection-1")
    assert radial.get_linewidth() > projection.get_linewidth()
    assert any("vᵣ = v∥" in text.get_text() for text in axis.texts)


def test_guided_arrows_are_drawn_on_the_gray_packet_ray() -> None:
    """The displayed arrows share the ray coordinates without screen offsets."""
    for case in (
        analyze_case(linked_inputs()),
        analyze_case(ExplorerInputs(1, 3, 1, 3, 0, -3, -2)),
    ):
        figure = rendering.render_guided_case(case)
        figure.canvas.draw()
        axis = next(
            item
            for item in figure.axes
            if item.get_gid() == "guided vector construction"
        )
        trajectory = artist(axis, "vector-active-trajectory")
        assert trajectory.get_color() != rendering.INK
        for number in range(1, len(case.crossings) + 1):
            for name in ("vector-radial", "vector-projection"):
                arrow = artist(axis, f"{name}-{number}")
                if not isinstance(arrow, FancyArrowPatch):
                    continue
                point, _ = arrow_endpoints(axis, f"{name}-{number}")
                np.testing.assert_allclose(
                    arrow.get_transform().transform(point),
                    axis.transData.transform(point),
                    rtol=0,
                    atol=1e-10,
                )
                assert trajectory.get_linewidth() < arrow.get_linewidth()


def test_endpoint_cues_are_gray_and_behind_velocity_arrows() -> None:
    """Start and end remain readable context when a resonance nears the end."""
    case = analyze_case(linked_inputs(x_0=1.2))
    guided = rendering.render_guided_case(case)
    guided_axis = next(
        axis
        for axis in guided.axes
        if axis.get_gid() == "guided vector construction"
    )
    for axis in (vector_axis(case), guided_axis):
        arrow = artist(axis, "vector-projection-1")
        for label, coordinate in (
            ("start", case.inputs.x_0),
            ("end", case.x_end),
        ):
            text = next(item for item in axis.texts if item.get_text() == label)
            np.testing.assert_allclose(
                to_rgba(text.get_color()), to_rgba(rendering.MUTED)
            )
            assert text.get_zorder() < arrow.get_zorder()
            marker = next(
                item
                for item in axis.collections
                if len(item.get_offsets()) == 1
                and np.allclose(
                    item.get_offsets()[0], [coordinate, case.inputs.b]
                )
            )
            np.testing.assert_allclose(
                marker.get_edgecolors()[0], to_rgba(rendering.MUTED)
            )
            assert marker.get_zorder() < arrow.get_zorder()


def test_tangent_slots_make_one_construction_per_distinct_resonance() -> None:
    """Claim: a double algebraic root is one packet encounter and one vector.

    Regime: inbound crest touch at x=-0.9 plus a separate outbound crossing.
    Verification: the analysis has two distinct crossings, so the diagram has
    exactly two crossing points regardless of the four production root slots.
    """
    case = analyze_case(linked_inputs(w=0.50625))
    assert len(case.crossings) == 2
    axis = vector_axis(case)
    points = [
        item
        for item in axis.get_children()
        if (item.get_gid() or "").startswith("vector-point-")
    ]
    assert len(points) == 2
    np.testing.assert_allclose(
        artist(axis, "vector-point-1").get_offsets()[0][0], -0.9, atol=1e-12
    )


def test_target_miss_and_earlier_stop_have_no_misleading_vectors() -> None:
    """Claim: only distinct verified resonances on the active flight get vectors.

    Regimes: target above active range and an electron scatter before the first
    linked crossing. Verification: no construction point is present, while the
    full geometric ray and active segment remain visible for context.
    """
    for inputs in (
        linked_inputs(w=0.7),
        linked_inputs(s_stop=0.5),
    ):
        case = analyze_case(inputs)
        assert case.crossings == () or all(
            crossing.s > case.active_end for crossing in case.crossings
        )
        axis = vector_axis(case)
        assert not any(
            (item.get_gid() or "").startswith("vector-point-")
            for item in axis.get_children()
        )
        artist(axis, "vector-ray")
        artist(axis, "vector-active-trajectory")
        assert not any(isinstance(patch, Circle) for patch in axis.patches)


def test_guided_diagram_fits_between_shell_and_radial_law() -> None:
    """The compact construction sits in the right column at readable size."""
    case = analyze_case(linked_inputs())
    figure = rendering.render_guided_case(case)
    figure.canvas.draw()
    axes = {axis.get_gid(): axis for axis in figure.axes}
    shell = axes["guided shell geometry"]
    vectors = axes["guided vector construction"]
    radial = axes["guided radial law"]
    assert radial.get_position().y1 < vectors.get_position().y0
    assert vectors.get_position().y1 < shell.get_position().y0
    np.testing.assert_allclose(
        vectors.get_position(original=True).x0,
        shell.get_position(original=True).x0,
    )
    np.testing.assert_allclose(
        vectors.get_position(original=True).width,
        radial.get_position(original=True).width,
    )
    assert 300 <= vectors.get_window_extent().width <= 400
    assert vectors.get_window_extent().height >= 200
    assert figure.get_figheight() <= 9


def test_nearby_root_numbers_do_not_overlap_in_vector_panel() -> None:
    """Two close simple roots retain independently readable number labels."""
    case = analyze_case(linked_inputs(w=0.50624))
    assert len(case.crossings) == 3
    axis = vector_axis(case)
    labels = [text for text in axis.texts if text.get_text() in {"1", "2", "3"}]
    assert len(labels) == 3
    renderer = axis.figure.canvas.get_renderer()
    extents = [label.get_window_extent(renderer) for label in labels]
    for first_idx, first in enumerate(extents):
        for second in extents[first_idx + 1 :]:
            assert not first.overlaps(second)


def test_detailed_figure_also_contains_shell_free_vector_construction() -> None:
    """The detailed export includes the new geometric explanation."""
    case = analyze_case(linked_inputs())
    figure = rendering.render_case(case)
    figure.canvas.draw()
    axes = {axis.get_gid(): axis for axis in figure.axes}
    vectors = axes["vector construction"]
    assert not any(isinstance(patch, Circle) for patch in vectors.patches)
    for number in (1, 2, 3):
        artist(vectors, f"vector-point-{number}")
