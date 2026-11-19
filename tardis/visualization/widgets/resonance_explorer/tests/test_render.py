"""Image regression checks for every lesson plot."""

from pathlib import Path

import pytest
from matplotlib.testing.compare import compare_images

from tardis.visualization.widgets.resonance_explorer.app import (
    ResonanceExplorerWidget,
)
from tardis.visualization.widgets.resonance_explorer.lessons import LESSONS


@pytest.mark.parametrize(
    ("lesson_key", "state_key"),
    [(lesson.key, state.key) for lesson in LESSONS for state in lesson.states],
)
def test_lesson_plot_matches_reference_image(
    lesson_key: str,
    state_key: str,
    regression_data: object,
    tmp_path: Path,
) -> None:
    """Claim: the plot a reader sees for each lesson step is unchanged."""
    widget = ResonanceExplorerWidget()
    widget.lesson_select.value = lesson_key
    widget.lesson_state_select.value = state_key
    figure = widget.figure_pane.object
    image_name = f"{regression_data.fname_prefix}.png"
    expected = regression_data.absolute_regression_data_dir / image_name

    if regression_data.enable_generate_reference:
        expected.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(expected)
        pytest.skip("Generated reference image")
    actual = tmp_path / image_name
    figure.savefig(actual)
    assert compare_images(str(expected), str(actual), tol=1e-3) is None
