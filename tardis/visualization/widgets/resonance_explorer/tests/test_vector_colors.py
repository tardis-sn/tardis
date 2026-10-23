"""Color consistency between the shell field and resonance radial vectors."""

import math
from collections.abc import Callable

import numpy as np
import pytest
from matplotlib.colors import TwoSlopeNorm, to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch
from matplotlib.patheffects import Stroke

from tardis.visualization.widgets.resonance_explorer import (
    CaseAnalysis,
    ExplorerInputs,
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer import render as rendering


@pytest.fixture
def linked_case() -> CaseAnalysis:
    """Give three resonances spanning blue and red shell-field regions."""
    inputs = ExplorerInputs(
        r_in=1.0,
        r_out=3.0,
        v_in=-1.34375,
        v_out=0.65625,
        b=1.2,
        x_0=-math.sqrt(7.56),
        w=0.3,
    )
    case = analyze_case(inputs)
    assert len(case.crossings) == 3
    return case


@pytest.mark.parametrize(
    ("renderer", "axis_gid"),
    [
        (rendering.render_guided_case, "guided vector construction"),
    ],
)
def test_radial_arrow_colors_follow_the_shell_field(
    linked_case: CaseAnalysis,
    renderer: Callable[..., Figure],
    axis_gid: str,
) -> None:
    """Claim: a resonance radial arrow has its local shell-field color.

    Regime: one oblique flight crosses both signs of the linear radial law.
    Verification: interpolate velocity directly between shell-edge inputs,
    then apply the field's endpoint-centered color scale at that radius.
    The edge interpolation is independent of the renderer's ``m*r + q``
    calculation. Colors are exact deterministic RGBA values.
    """
    figure = renderer(linked_case)
    axis = next(ax for ax in figure.axes if ax.get_gid() == axis_gid)
    inputs = linked_case.inputs
    velocity_scale = max(abs(inputs.v_in), abs(inputs.v_out))
    color_norm = TwoSlopeNorm(
        vmin=-velocity_scale, vcenter=0.0, vmax=velocity_scale
    )
    for number, crossing in enumerate(linked_case.crossings, start=1):
        assert crossing.x is not None
        radius = math.hypot(crossing.x, inputs.b)
        fraction = (radius - inputs.r_in) / (inputs.r_out - inputs.r_in)
        radial_velocity = inputs.v_in + fraction * (inputs.v_out - inputs.v_in)
        expected_color = rendering.PALETTE(color_norm(radial_velocity))
        arrow = next(
            artist
            for artist in axis.get_children()
            if artist.get_gid() == f"vector-radial-{number}"
        )
        assert isinstance(arrow, FancyArrowPatch)
        np.testing.assert_allclose(
            arrow.get_facecolor(), expected_color, rtol=0, atol=1e-12
        )


def test_radial_arrows_have_a_very_thin_black_outline(
    linked_case: CaseAnalysis,
) -> None:
    """Each field-colored radial arrow has a visible black path outline."""
    figure = rendering.render_guided_case(linked_case)
    axis = next(
        ax for ax in figure.axes if ax.get_gid() == "guided vector construction"
    )
    for number in range(1, len(linked_case.crossings) + 1):
        arrow = next(
            artist
            for artist in axis.get_children()
            if artist.get_gid() == f"vector-radial-{number}"
        )
        assert isinstance(arrow, FancyArrowPatch)
        strokes = [
            effect
            for effect in arrow.get_path_effects()
            if isinstance(effect, Stroke)
        ]
        assert strokes, f"radial arrow {number} has no black outline"
        assert any(
            np.allclose(to_rgba(stroke._gc["foreground"]), to_rgba("black"))
            and 0 < stroke._gc["linewidth"] - arrow.get_linewidth() <= 1.5
            for stroke in strokes
        ), f"radial arrow {number} needs a very thin black outline"


def test_zero_radial_speed_marker_uses_the_midpoint_color(
    linked_case: CaseAnalysis,
) -> None:
    """A stationary radial velocity keeps the field's near-white zero color.

    Regime: line target w=0 gives exact stagnation resonances on both sides
    of closest approach. Verification: the zero-vector glyph has the same
    colormap midpoint as the shell's stagnation ring, with a dark outline so
    this pale color remains visible.
    """
    inputs = linked_case.inputs
    case = analyze_case(
        ExplorerInputs(
            r_in=inputs.r_in,
            r_out=inputs.r_out,
            v_in=inputs.v_in,
            v_out=inputs.v_out,
            b=inputs.b,
            x_0=inputs.x_0,
            w=0.0,
        )
    )
    figure = rendering.render_guided_case(case)
    axis = next(
        ax for ax in figure.axes if ax.get_gid() == "guided vector construction"
    )
    zero_color = rendering.PALETTE(0.5)
    for number in (1, 3):
        crossing = case.crossings[number - 1]
        assert crossing.x is not None
        radius = math.hypot(crossing.x, inputs.b)
        np.testing.assert_allclose(
            inputs.v_in
            + (radius - inputs.r_in)
            * (inputs.v_out - inputs.v_in)
            / (inputs.r_out - inputs.r_in),
            0.0,
            atol=1e-12,
        )
        marker = next(
            artist
            for artist in axis.get_children()
            if artist.get_gid() == f"vector-radial-{number}"
        )
        colors = [*marker.get_facecolors(), *marker.get_edgecolors()]
        assert any(np.allclose(color, zero_color) for color in colors)
        has_black_edge = any(
            np.allclose(color, to_rgba("black"))
            for color in marker.get_edgecolors()
        )
        has_black_stroke = any(
            isinstance(effect, Stroke)
            and np.allclose(to_rgba(effect._gc["foreground"]), to_rgba("black"))
            for effect in marker.get_path_effects() or ()
        )
        assert has_black_edge or has_black_stroke
