"""Headless checks of the scientific coordinates shown by the linked figure."""

import math

import numpy as np
from matplotlib.colors import to_rgba

from tardis.visualization.widgets.resonance_explorer import (
    ExplorerInputs,
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.render import (
    render_case,
    render_guided_case,
)


def linked_inputs(**changes: float) -> ExplorerInputs:
    """Return the three-crossing flyby in illustrative units."""
    values = dict(
        r_in=1.0,
        r_out=3.0,
        v_in=-1.34375,
        v_out=0.65625,
        b=1.2,
        x_0=-math.sqrt(3**2 - 1.2**2),
        w=0.3,
    )
    values.update(changes)
    return ExplorerInputs(**values)


def find_axis(figure: object, name: str) -> object:
    """Find one named view without relying on the colorbar position."""
    return next(ax for ax in figure.axes if ax.get_gid() == name)


def find_artist(axis: object, gid: str) -> object:
    """Find a scientific event marker or active segment by its stable id."""
    return next(
        artist for artist in axis.get_children() if artist.get_gid() == gid
    )


def test_linked_crossings_align_with_geometry_and_target() -> None:
    """Claim: the figure links published roots across physical and curve views.

    Regime: complete b=1.2 inner-miss flyby, m=1, q=-2.34375, w=0.3.
    Verification: published unsquared roots fix x; geometry has y=b, and
    the projected curve meets the independently specified horizontal w.
    The 1e-6 x tolerance reflects the six printed decimals of those roots.
    """
    case = analyze_case(linked_inputs())
    figure = render_case(case)
    figure.canvas.draw()
    geometry = find_axis(figure, "shell geometry")
    projected = find_axis(figure, "projected velocity")
    gradient = find_axis(figure, "projected gradient")
    for candidate in case.candidates:
        if candidate.x is not None:
            raw_slot = find_artist(projected, f"raw-slot-{candidate.slot_idx}")
            np.testing.assert_allclose(
                raw_slot.get_offsets()[0][0], candidate.x
            )
    for crossing_idx, expected_x in enumerate(
        [-1.555887, -0.341473, 2.395531], start=1
    ):
        marker_id = f"crossing-{crossing_idx}"
        geometry_position = find_artist(geometry, marker_id).get_offsets()[0]
        projected_position = find_artist(projected, marker_id).get_offsets()[0]
        np.testing.assert_allclose(
            geometry_position, [expected_x, 1.2], atol=1e-6
        )
        np.testing.assert_allclose(
            projected_position, [expected_x, 0.3], atol=1e-6
        )
        geometry_pixel_x = geometry.transData.transform(geometry_position)[0]
        projected_pixel_x = projected.transData.transform(projected_position)[0]
        gradient_pixel_x = gradient.transData.transform((expected_x, 0))[0]
        np.testing.assert_allclose(
            [geometry_pixel_x, projected_pixel_x, gradient_pixel_x],
            geometry_pixel_x,
            atol=0.5,
        )
    assert geometry.get_ylim()[0] < -3
    assert geometry.get_ylim()[1] > 3


def test_target_miss_has_no_crossing_markers() -> None:
    """Claim: a target above the active projected range shows no fake roots.

    Verification: w=0.7 exceeds the linked crest near 0.50625.
    This is a rendering check; the analyzer determines the root count.
    """
    case = analyze_case(linked_inputs(w=0.7))
    figure = render_case(case)
    figure.canvas.draw()
    for name in ("shell geometry", "projected velocity"):
        axis = find_axis(figure, name)
        assert not any(
            artist.get_gid() and artist.get_gid().startswith("crossing-")
            for artist in axis.get_children()
        )
    target_axis = find_axis(figure, "projected velocity")
    assert "w = 0.7" in target_axis.get_legend().get_texts()[0].get_text()


def test_unverified_production_slot_has_distinct_plot_disposition() -> None:
    """A production-accepted root with a large unsquared error is visible but unnumbered."""
    case = analyze_case(
        ExplorerInputs(
            1,
            3,
            -1103.1408330241957,
            -1103.137596530609,
            1.0987592278970888,
            -2.7915458368278663,
            22.32587042642611,
        )
    )
    figure = render_case(case)
    figure.canvas.draw()
    projected = find_axis(figure, "projected velocity")
    slot = find_artist(projected, "candidate-slot-1")
    np.testing.assert_allclose(slot.get_edgecolors()[0], to_rgba("#a44d77"))
    assert not any(
        artist.get_gid() and artist.get_gid().startswith("crossing-")
        for artist in projected.get_children()
    )
    assert "production accepted, unsquared unverified" in (
        label.get_text()
        for legend in figure.legends
        for label in legend.get_texts()
    )


def test_nearest_verified_line_has_its_own_plot_marker() -> None:
    """The next eligible line is identifiable without reading the status text."""
    case = analyze_case(linked_inputs())
    figure = render_case(case)
    figure.canvas.draw()
    for name in ("shell geometry", "projected velocity"):
        axis = find_axis(figure, name)
        marker = find_artist(axis, "nearest-verified-line")
        np.testing.assert_allclose(
            marker.get_offsets()[0][0], case.nearest_resonance.x
        )
    assert "nearest verified line" in (
        label.get_text()
        for legend in figure.legends
        for label in legend.get_texts()
    )


def test_close_crossing_number_labels_do_not_overlap() -> None:
    """Two nearby eligible roots retain independently legible numbers."""
    case = analyze_case(linked_inputs(w=0.50624))
    assert len(case.crossings) == 3
    figure = render_case(case)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for name in ("shell geometry", "projected velocity"):
        axis = find_axis(figure, name)
        labels = [
            label for label in axis.texts if label.get_text() in {"1", "2", "3"}
        ]
        assert len(labels) == 3
        extents = [label.get_window_extent(renderer) for label in labels]
        assert not extents[0].overlaps(extents[1])


def test_geometric_only_turn_has_distinct_plot_marker() -> None:
    """A tolerance-excluded turn remains visible without implying a packet split."""
    impact = 1.2
    start = -math.sqrt(9 - impact**2)
    turn_x = start + 5e-15
    intercept = -((turn_x**2 + impact**2) ** 1.5) / impact**2
    case = analyze_case(
        ExplorerInputs(1, 3, 1 + intercept, 3 + intercept, impact, start, 0)
    )
    assert case.turns and not case.transport_turns
    figure = render_case(case)
    figure.canvas.draw()
    geometry = find_axis(figure, "shell geometry")
    for turn_idx, turn in enumerate(case.turn_candidates):
        if turn.active and not turn.transport_active:
            marker = find_artist(geometry, f"turn-{turn_idx}")
            np.testing.assert_allclose(
                marker.get_edgecolors()[0], to_rgba("#a66b2e")
            )
    assert "geometric turn, transport excluded" in (
        label.get_text()
        for legend in figure.legends
        for label in legend.get_texts()
    )


def test_earlier_stop_shortens_only_the_active_segment() -> None:
    """Claim: an assumed scatter ends the bold packet segment at x0+s_stop.

    Regime: linked ray stopped two illustrative length units after launch.
    Verification: the input path equation x=x0+s fixes the endpoint.
    """
    inputs = linked_inputs(s_stop=2.0)
    figure = render_case(analyze_case(inputs))
    figure.canvas.draw()
    active_ray = find_artist(find_axis(figure, "shell geometry"), "active-ray")
    active_projection = find_artist(
        find_axis(figure, "projected velocity"), "active-projection"
    )
    np.testing.assert_allclose(
        active_ray.get_xdata(), [inputs.x_0, inputs.x_0 + 2]
    )
    np.testing.assert_allclose(
        [active_projection.get_xdata()[0], active_projection.get_xdata()[-1]],
        [inputs.x_0, inputs.x_0 + 2],
    )


def test_exact_turn_target_gets_contact_marker() -> None:
    """Claim: the w=0.50625 crest contact is visible as a target touch.

    Verification: analytic G=0 at x=-0.9 and v_parallel(-0.9)=0.50625.
    The marker is separate from solver acceptance of a repeated root.
    """
    figure = render_case(analyze_case(linked_inputs(w=0.50625)))
    figure.canvas.draw()
    projected = find_axis(figure, "projected velocity")
    contacts = [
        artist
        for artist in projected.get_children()
        if artist.get_gid() and artist.get_gid().startswith("target-contact-")
    ]
    assert len(contacts) == 1
    np.testing.assert_allclose(contacts[0].get_offsets()[0], [-0.9, 0.50625])


def test_radial_inward_active_curve_avoids_origin_singularity() -> None:
    """Claim: a radial inward shell segment remains finite in all linked views.

    Regime: b=0 and x from -3 to the first inner hit at -1.
    Verification: the flight never enters the excised origin, even though
    the complete reference chord passes through it.
    """
    case = analyze_case(ExplorerInputs(1, 3, -1.34375, 0.65625, 0, -3, 0.3))
    figure = render_case(case)
    figure.canvas.draw()
    for axis_name, gid in (
        ("projected velocity", "active-projection"),
        ("projected gradient", "active-gradient"),
    ):
        curve = find_artist(find_axis(figure, axis_name), gid)
        assert np.all(np.isfinite(curve.get_ydata()))


def test_selected_raw_slot_is_highlighted_in_all_linked_views() -> None:
    """Claim: inspecting a finite solver slot rings that slot's same x.

    This is a UI wiring check: selection must not change the analysis.
    """
    case = analyze_case(linked_inputs())
    selected = next(
        candidate for candidate in case.candidates if candidate.x is not None
    )
    figure = render_case(case, selected_slot_index=selected.slot_idx)
    figure.canvas.draw()
    for name in ("shell geometry", "projected velocity", "projected gradient"):
        axis = find_axis(figure, name)
        ring = find_artist(axis, f"selected-slot-{selected.slot_idx}")
        np.testing.assert_allclose(ring.get_offsets()[0][0], selected.x)


def test_guided_homology_comparison_keeps_axes_and_target_fixed() -> None:
    """Changing only q preserves geometry, target, and projected plot scale."""
    linked = analyze_case(linked_inputs())
    homology = analyze_case(linked_inputs(v_in=1, v_out=3))
    first = render_guided_case(homology, reference=linked)
    second = render_guided_case(linked, reference=homology)
    for name in (
        "guided projected velocity",
        "guided shell geometry",
        "guided radial law",
        "guided projected gradient",
    ):
        first_axis = find_axis(first, name)
        second_axis = find_axis(second, name)
        np.testing.assert_allclose(
            first_axis.get_xlim(), second_axis.get_xlim()
        )
        np.testing.assert_allclose(
            first_axis.get_ylim(), second_axis.get_ylim()
        )
    projection = find_axis(second, "guided projected velocity")
    gradient = find_axis(second, "guided projected gradient")
    assert "G>0" in gradient.get_title(loc="left")
    assert "G<0" in gradient.get_title(loc="left")
    assert "homology reference" in {
        label.get_text() for label in projection.get_legend().get_texts()
    }
    assert any(
        artist.get_gid() == "active-flight-shading"
        for artist in projection.patches
    )
    assert any(
        line.get_gid() == "line-target"
        and np.all(np.asarray(line.get_ydata()) == 0.3)
        for line in projection.lines
    )
    assert len(linked.crossings) == 3
    for number, candidate in enumerate(linked.crossings, start=1):
        point = find_artist(projection, f"crossing-{number}")
        np.testing.assert_allclose(point.get_offsets()[0], [candidate.x, 0.3])


def test_guided_selected_slot_and_collapse_crossings_are_visible() -> None:
    """The compact view identifies two inbound inflow roots and a raw slot."""
    case = analyze_case(linked_inputs(v_in=-2.5, v_out=-0.5, w=0.6))
    assert len(case.crossings) == 2
    selected = next(c for c in case.candidates if c.x is not None)
    figure = render_guided_case(case, selected_slot_index=selected.slot_idx)
    projected = find_axis(figure, "guided projected velocity")
    radial = find_axis(figure, "guided radial law")
    assert max(radial.lines[1].get_ydata()) < 0
    find_artist(projected, f"selected-slot-{selected.slot_idx}")
    for number in (1, 2):
        find_artist(projected, f"crossing-{number}")
