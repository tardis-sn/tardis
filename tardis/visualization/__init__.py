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
from tardis.visualization.sankey import (
    plot_packet_history_sankey,
    plot_multi_packet_sankey_after_carbon,
)


print("Initializing tabulator and plotly panel extensions for widgets to work")
import panel as pn
pn.extension("tabulator", "plotly")