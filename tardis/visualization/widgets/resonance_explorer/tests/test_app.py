"""Behavior checks for the explorer interface."""

from bokeh.document import Document

from tardis.visualization.widgets.resonance_explorer.app import (
    ResonanceExplorerWidget,
)


def test_default_case_serializes_to_a_browser_document() -> None:
    """Claim: the layout contains only values Panel can serialize."""
    widget = ResonanceExplorerWidget()
    document = Document()

    widget.layout.server_doc(document)

    assert document.roots


def test_line_velocity_slider_refills_the_same_plot() -> None:
    """Claim: a slider step sends new data into the existing plot."""
    widget = ResonanceExplorerWidget()
    plot = widget.figure_pane.object

    widget.sliders["line_velocity"].value = 0.5

    assert widget.analysis.inputs.line_velocity == 0.5
    assert widget.figure_pane.object is plot
    assert plot.select_one({"name": "line-target"}).data["y"] == [0.5]
