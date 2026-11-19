"""The explorer's notebook entry point follows other Panel widgets."""

import panel as pn

from tardis.visualization import ResonanceExplorerWidget


def test_explorer_widget_displays_a_live_panel_layout() -> None:
    """A user can construct and display the explorer from tardis.visualization."""
    widget = ResonanceExplorerWidget()
    displayed = widget.display()

    assert isinstance(displayed, pn.layout.Row)
    assert displayed is widget.layout
    assert widget.analysis is not None
