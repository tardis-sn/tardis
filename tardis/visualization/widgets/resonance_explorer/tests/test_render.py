"""Regression checks for the plotted data of every lesson view."""

import json

import numpy as np
import pytest
from bokeh.models import ColumnDataSource, Label, Range1d

from tardis.conftest import sync_ndarray_assert_allclose
from tardis.visualization.widgets.resonance_explorer.app import (
    ResonanceExplorerWidget,
)
from tardis.visualization.widgets.resonance_explorer.lessons import LESSONS


@pytest.mark.parametrize(
    ("lesson_key", "state_key"),
    [(lesson.key, state.key) for lesson in LESSONS for state in lesson.states],
)
def test_lesson_plot_matches_regression_data(
    lesson_key: str, state_key: str, regression_data: object
) -> None:
    """Claim: the curves, markers, and labels a reader sees are unchanged."""
    widget = ResonanceExplorerWidget()
    widget.lesson_select.value = lesson_key
    widget.lesson_state_select.value = state_key
    plot = widget.figure_pane.object
    numbers = []
    styles = {}
    # Unnamed sources only hold the empty legend key lines
    named_sources = [
        source
        for source in plot.select({"type": ColumnDataSource})
        if source.name
    ]
    for source in sorted(named_sources, key=lambda item: item.name):
        for column, values in sorted(source.data.items()):
            array = np.asarray(values)
            if array.dtype.kind in "fiu":
                numbers.append(array)
            else:
                styles[f"{source.name}.{column}"] = list(values)
    numbers.extend(
        [plot_range.start, plot_range.end]
        for plot_range in sorted(
            plot.select({"type": Range1d}), key=lambda item: item.name
        )
    )
    for label in plot.select({"type": Label}):
        styles[label.name] = [label.text, label.visible]

    sync_ndarray_assert_allclose(regression_data, *numbers)
    text = json.dumps(styles, indent=1, sort_keys=True, ensure_ascii=False)
    assert regression_data.sync_str(text) == text
