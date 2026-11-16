"""The homology comparison changes overlays without rescaling a guided view."""

import numpy as np
import numpy.testing as npt
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch

from tardis.visualization.widgets.resonance_explorer.app import (
    ResonanceExplorerWidget,
)


def _plot_scales_and_shell_colors(
    figure: Figure,
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], np.ndarray, np.ndarray]:
    """Capture linked axes, shell colors, and radial-arrow colors."""
    names = (
        "guided projected velocity",
        "guided projected gradient",
        "guided radial law",
        "guided shell geometry",
    )
    scales = {
        name: (
            np.asarray(
                next(
                    ax for ax in figure.axes if ax.get_gid() == name
                ).get_xlim()
            ),
            np.asarray(
                next(
                    ax for ax in figure.axes if ax.get_gid() == name
                ).get_ylim()
            ),
        )
        for name in names
    }
    shell_axis = next(
        ax for ax in figure.axes if ax.get_gid() == "guided shell geometry"
    )
    shell_colors = np.asarray(
        [
            patch.get_facecolor()
            for patch in shell_axis.patches
            if isinstance(patch, Circle)
        ]
    )
    vector_axis = next(
        ax for ax in figure.axes if ax.get_gid() == "guided vector construction"
    )
    arrow_colors = np.asarray(
        [
            patch.get_facecolor()
            for patch in vector_axis.patches
            if isinstance(patch, FancyArrowPatch)
            and (patch.get_gid() or "").startswith("vector-radial-")
        ]
    )
    return scales, shell_colors, arrow_colors


@pytest.mark.parametrize("lesson_key", ["homology", "linked"])
def test_comparison_toggle_preserves_guided_scales_and_shell_colors(
    lesson_key: str,
) -> None:
    """Only the dashed reference changes when the paired overlay is toggled."""
    application = ResonanceExplorerWidget()
    application.lesson_select.value = lesson_key
    application.comparison_toggle.value = False
    hidden = application.figure_pane.object
    assert isinstance(hidden, Figure)
    hidden_scales, hidden_colors, hidden_arrow_colors = (
        _plot_scales_and_shell_colors(hidden)
    )
    assert not any(
        "reference" in line.get_label()
        for ax in hidden.axes
        for line in ax.lines
    )

    application.comparison_toggle.value = True
    shown = application.figure_pane.object
    assert isinstance(shown, Figure)
    shown_scales, shown_colors, shown_arrow_colors = (
        _plot_scales_and_shell_colors(shown)
    )
    assert any(
        "reference" in line.get_label()
        for ax in shown.axes
        for line in ax.lines
    )

    for name in hidden_scales:
        for before, after in zip(
            hidden_scales[name], shown_scales[name], strict=True
        ):
            npt.assert_array_equal(
                before, after, err_msg=f"{lesson_key}: {name}"
            )
    npt.assert_array_equal(
        hidden_colors,
        shown_colors,
        err_msg=f"{lesson_key}: shell colormap must retain one normalization",
    )
    npt.assert_array_equal(
        hidden_arrow_colors,
        shown_arrow_colors,
        err_msg=f"{lesson_key}: radial arrows must retain shell colors",
    )


def test_paired_guided_lessons_share_axes_without_reference_overlay() -> None:
    """The reference and linked cases retain their scales when switched."""
    application = ResonanceExplorerWidget()
    application.comparison_toggle.value = False
    reference = application.figure_pane.object
    assert isinstance(reference, Figure)
    reference_scales, _, _ = _plot_scales_and_shell_colors(reference)

    application.lesson_select.value = "linked"
    linked = application.figure_pane.object
    assert isinstance(linked, Figure)
    linked_scales, _, _ = _plot_scales_and_shell_colors(linked)
    for name in reference_scales:
        for before, after in zip(
            reference_scales[name], linked_scales[name], strict=True
        ):
            npt.assert_array_equal(before, after, err_msg=name)
