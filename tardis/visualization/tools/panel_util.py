"""Panel utility functions for visualization tools."""

# Panel extension initialization - lazily initialized to avoid output at import time
_panel_initialized = False


def ensure_panel_extension():
    """
    Lazily initialize Panel extensions for tabulator and plotly.

    This function should be called before using Panel widgets that require
    these extensions. It only initializes once per session.
    """
    global _panel_initialized
    if not _panel_initialized:
        import panel as pn

        pn.extension("tabulator", "plotly")
        _panel_initialized = True
