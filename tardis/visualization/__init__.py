"""Visualization tools and widgets for TARDIS."""

import logging

import panel as pn

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
    "LIVPlotter",
    "RPacketPlotter",
    "SDECPlotter",
    "CustomAbundanceWidget",
    "GrotrianWidget",
    "LineInfoWidget",
    "shell_info_from_hdf",
    "shell_info_from_simulation",
]

logger = logging.getLogger(__name__)
logger.info("Initializing tabulator and plotly panel extensions for widgets to work")

pn.extension("tabulator", "plotly")