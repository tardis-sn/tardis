"""Headless checks of the scientific coordinates shown by the linked figure."""

import math

import numpy as np

from tardis.visualization.widgets.resonance_explorer import (
    ExplorerInputs,
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.render import (
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
