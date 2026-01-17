"""Visualization tools and widgets for TARDIS."""

from tardis.visualization.tools.convergence_plot import ConvergencePlots
from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.visualization.tools.rpacket_plot import RPacketPlotter
from tardis.visualization.tools.sdec_plot import SDECPlotter
from tardis.visualization.widgets.custom_abundance import CustomAbundanceWidget
from tardis.visualization.widgets.grotrian import GrotrianWidget
from tardis.visualization.widgets.line_info import LineInfoWidget
from tardis.visualization.widgets.shell_info import (
    shell_info_from_hdf,
    shell_info_from_simulation,
)



__all__ = [
    "ConvergencePlots",
    "CustomAbundanceWidget",
    "GrotrianWidget",
    "LIVPlotter",
    "LineInfoWidget",
    "RPacketPlotter",
    "SDECPlotter",
    "shell_info_from_hdf",
    "shell_info_from_simulation",
]

_panel_initialized = False


def _ensure_panel_initialized():
    """initialize panel extensions lazily to prevent import side effects."""
    global _panel_initialized
    if not _panel_initialized:
        import panel as pn
        pn.extension("tabulator", "plotly")
        _panel_initialized = True
