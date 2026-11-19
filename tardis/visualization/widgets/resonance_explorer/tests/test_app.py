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


def test_line_velocity_slider_redraws_plot() -> None:
    """Claim: moving a slider updates the inputs and replaces the plot."""
    widget = ResonanceExplorerWidget()
    previous_figure = widget.figure_pane.object

    widget.sliders["w"].value = 0.5

    assert widget.analysis.inputs.w == 0.5
    assert widget.figure_pane.object is not previous_figure
