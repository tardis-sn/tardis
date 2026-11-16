"""Visualization tools and widgets for TARDIS."""

import panel as pn

from tardis.visualization.tools.convergence_plot import ConvergencePlots
from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.visualization.tools.rpacket_plot import RPacketPlotter
from tardis.visualization.tools.sdec_plot import SDECPlotter
from tardis.visualization.widgets.custom_abundance import CustomAbundanceWidget
from tardis.visualization.widgets.grotrian import GrotrianWidget
from tardis.visualization.widgets.line_info import LineInfoWidget
from tardis.visualization.widgets.resonance_explorer import (
    ResonanceExplorerWidget,
)
from tardis.visualization.widgets.shell_info import ShellInfoWidget

__all__ = [
    "ConvergencePlots",
    "CustomAbundanceWidget",
    "GrotrianWidget",
    "LIVPlotter",
    "LineInfoWidget",
    "RPacketPlotter",
    "ResonanceExplorerWidget",
    "SDECPlotter",
    "ShellInfoWidget",
]

print("Initializing tabulator and plotly panel extensions for widgets to work")
pn.extension("tabulator", "plotly", "katex")
