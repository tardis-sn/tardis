from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
from IPython.display import display
import matplotlib as mpl
import ipywidgets as widgets
from contextlib import suppress
from traitlets import TraitError


def transition_colors(iterations, name="jet"):
    """
    Function to create colorscale for convergence plots, returns a list of colors.

    Parameters
    ----------
    iterations : int
        Number of iterations.
    name : string
        Name of the colorscale. Defaults to "jet".

    Returns
    -------
    colors: list
    """
    cmap = mpl.cm.get_cmap(name, iterations)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(mpl.colors.rgb2hex(rgb))
    return colors


class ConvergencePlots(object):
    """
    Class to create and update convergence plots for visualizing convergence of the
    simulation.

    Parameters
    ----------
    iteration : int
        iteration number
    """

    def __init__(self, iterations, **kwargs):
        self.iterable_data = {}
        self.value_data = defaultdict(list)
        self.iterations = iterations
        self.current_iteration = 1
        self.luminosities = ["Emitted", "Absorbed", "Requested"]
        self.plasma_plot = None
        self.luminosity_plot = None

        if "colorscale" in kwargs:
            self.colorscale = transition_colors(
                kwargs["colorscale"], iterations=self.iterations
            )
        else:
            self.colorscale = transition_colors(iterations=self.iterations)

        if "plasma_plot_config" in kwargs:
            if not isinstance(kwargs["plasma_plot_config"], dict):
                raise TypeError("Expected dict in plasma_plot_config argument")
            self.plasma_plot_config = kwargs["plasma_plot_config"]

        if "luminosity_plot_config" in kwargs:
            if not isinstance(kwargs["luminosity_plot_config"], dict):
                raise TypeError(
                    "Expected dict in luminosity_plot_config argument"
                )
            self.luminosity_plot_config = kwargs["luminosity_plot_config"]

    def fetch_data(self, name=None, value=None, item_type=None):
        """
        This allows user to fetch data from the Simulation class.
        This data is stored and used when an iteration is completed.

        Parameters
        ----------
        name : string
            name of the data
        value : string or array
            string or array of quantities
        item_type : string
            either iterable or value

        """
        if item_type == "iterable":
            self.iterable_data[name] = value
        if item_type == "value":
            self.value_data[name].append(value)

    def create_plasma_plot(self):
        """
        Creates an empty plasma plot.
        The default layout can be overridden by passing plasma_plot_config dictionary in the run_tardis function.
        """
        fig = go.FigureWidget().set_subplots(rows=1, cols=2, shared_xaxes=True)
        fig.add_scatter(row=1, col=1)
        fig.add_scatter(row=1, col=2)

        fig = fig.update_layout(
            xaxis={
                "tickformat": "g",
                "title": r"$Shell~~Velocity$",
            },
            xaxis2={
                "tickformat": "g",
                "title": r"$Shell~~Velocity$",
                "matches": "x",
            },
            yaxis={
                "tickformat": "g",
                "title": r"$T_{rad}\ [K]$",
                "tickmode": "auto",
                "nticks": 15,
            },
            yaxis2={
                "tickformat": "g",
                "title": r"$W$",
                "tickmode": "auto",
                "nticks": 15,
            },
            height=450,
            legend_title_text="Iterations",
            margin=dict(l=10, r=135, b=25, t=25, pad=0),
        )

        # allows overriding default layout
        if hasattr(self, "plasma_plot_config"):
            self.override_plot_parameters(fig, self.plasma_plot_config)
        self.plasma_plot = fig

    def create_luminosity_plot(self):
        """
        Creates an empty luminosity plot.
        The default layout can be overridden by passing luminosity_plot_config dictionary in the run_tardis function.
        """
        line_colors = ["orangered", "lightseagreen", "indigo"]

        fig = go.FigureWidget().set_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            row_heights=[0.25, 0.5, 0.25],
        )

        for luminosity, line_color in zip(self.luminosities, line_colors):
            fig.add_scatter(
                name=luminosity + "<br>Luminosity",
                mode="lines",
                row=2,
                col=1,
                marker_color=line_color,
                legendgroup=luminosity,
            )

        fig.add_scatter(
            name="Residual<br>Luminosity",
            row=3,
            col=1,
            marker_color="cornflowerblue",
            mode="lines",
        )

        fig.add_scatter(
            name="Inner<br>Boundary<br>Temperature",
            row=1,
            col=1,
            hovertext="text",
            marker_color="crimson",
            mode="lines",
        )

        # 3 y axes and 3 x axes correspond to the 3 subplots in the luminosity convergence plot
        fig = fig.update_layout(
            xaxis=dict(range=[0, self.iterations + 1], dtick=2),
            xaxis2=dict(
                matches="x",
                range=[0, self.iterations + 1],
                dtick=2,
            ),
            xaxis3=dict(
                matches="x2",
                title=r"$\mbox{Iteration Number}$",
                dtick=2,
            ),
            yaxis=dict(
                title=r"$\mbox{T}_{inner}\ [K]$",
                automargin=True,
                side="top",
                tickformat="g",
                exponentformat="e",
                tickmode="auto",
                nticks=4,
            ),
            yaxis2=dict(
                exponentformat="e",
                title=r"$\mbox{Luminosity}~(erg~s^{-1})$",
                title_font_size=13,
                automargin=True,
                tickmode="auto",
                nticks=7,
            ),
            yaxis3=dict(
                exponentformat="e",
                title=r"$~~\mbox{Residual}\\\mbox{Luminosity(%)}$",
                title_font_size=12,
                automargin=True,
                tickformat="%",
                tickmode="auto",
                nticks=4,
            ),
            legend_tracegroupgap=0,
            height=630,
            hoverlabel_align="right",
            legend_title_text="Luminosity",
            hoverlabel_font_color="white",
            margin=dict(b=25, t=25, pad=0),
        )

        # allows overriding default layout
        if hasattr(self, "luminosity_plot_config"):
            self.override_plot_parameters(fig, self.luminosity_plot_config)

        self.luminosity_plot = fig

    def override_plot_parameters(self, fig, parameters):
        """
        Overrides default plot properties.

        Parameters
        ----------
        fig : go.FigureWidget
            FigureWidget object to be updated
        parameters : dict
            Dictionary used to update the default plot style. The dictionary should have a structure
            like that of go.FigureWidget.to_dict(), for more information please see https://plotly.com/python/figure-structure/
        """
        for key, value in parameters.items():
            if key == "data":
                for trace in list(fig.data):
                    self.override_plot_parameters(trace, value)
            else:
                if type(value) == dict:
                    self.override_plot_parameters(fig[key], value)
                else:
                    fig[key] = value

    def build(self, display_plot=True):
        """
        Creates empty plasma and luminosity convergence plots and displays them together.

        Parameters
        ----------
        display_plot : bool
            Displays empty plots. Defaults to True.
        """
        self.create_plasma_plot()
        self.create_luminosity_plot()
        if display_plot:
            display(
                widgets.VBox(
                    [self.plasma_plot, self.luminosity_plot],
                )
            )

    def update_plasma_plots(self):
        """
        Updates plasma convergence plots every iteration.
        """
        x = self.iterable_data["velocity"].value.tolist()
        x = [item / 100000 for item in x]

        customdata = len(x) * [
            "<br>"
            + "Emitted Luminosity: "
            + f'{self.value_data["Absorbed"][-1]:.2g}'
            + "<br>"
            + "Requested Luminosity: "
            + f'{self.value_data["Requested"][-1]:.2g}'
            + "<br>"
            + "Absorbed Luminosity: "
            + f'{self.value_data["Requested"][-1]:.2g}'
        ]

        # this adds radiation temperature vs shell velocity plot
        self.plasma_plot.add_scatter(
            x=x,
            y=self.iterable_data["t_rad"],
            line_color=self.colorscale[self.current_iteration - 1],
            row=1,
            col=1,
            name=self.current_iteration,
            legendgroup=f"group-{self.current_iteration}",
            showlegend=False,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.2f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

        # this adds dilution factor vs shell velocity plot
        self.plasma_plot.add_scatter(
            x=x,
            y=self.iterable_data["w"],
            line_color=self.colorscale[self.current_iteration - 1],
            row=1,
            col=2,
            legendgroup=f"group-{self.current_iteration}",
            name=self.current_iteration,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.2f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

    def update_luminosity_plot(self):
        """
        Updates luminosity convergence plots every iteration.
        """
        x = list(range(1, self.iterations + 1))
        with self.luminosity_plot.batch_update():
            for index, luminosity in zip(range(3), self.luminosities):
                self.luminosity_plot.data[index].x = x
                self.luminosity_plot.data[index].y = self.value_data[luminosity]
                self.luminosity_plot.data[index].hovertemplate = (
                    "<b>%{y:.2g}</b>" + "<br>at X = %{x}<br>"
                )

                y = [
                    ((emitted - requested) * 100) / requested
                    for emitted, requested in zip(
                        self.value_data["Emitted"], self.value_data["Requested"]
                    )
                ]
                self.luminosity_plot.data[-2].x = x
                self.luminosity_plot.data[-2].y = y
                self.luminosity_plot.data[
                    -2
                ].hovertemplate = "Residual Luminosity: %{y:.2f}% at X = %{x:,.0f}<extra></extra>"

                self.luminosity_plot.data[-1].x = x
                self.luminosity_plot.data[-1].y = self.value_data["t_inner"]
                self.luminosity_plot.data[
                    -1
                ].hovertemplate = "Inner Body Temperature: %{y:.2f} at X = %{x:,.0f}<extra></extra>"

    def update(self, export_cplots=False, last=False):
        """
        Updates the plasma and the luminosity convergence plots every iteration
        and displays them again after the last iteration for exporting if desired.

        Parameters
        ----------
        export_cplots : bool
            Displays  the convergence plots using plotly's 'notebook_connected' renderer.
            This displays the plot in notebooks when shared on platforms like nbviewer.
            Please see https://plotly.com/python/renderers/ for more information.
            Defaults to False.
        last : bool
            True if it's last iteration.
        """
        if self.iterable_data != {}:
            if self.current_iteration == 1:
                self.build()

            self.update_plasma_plots()
            self.update_luminosity_plot()

            if hasattr(self, "plasma_plot_config") and last:
                if "data" in self.plasma_plot_config:
                    self.override_plot_parameters(
                        self.plasma_plot, self.plasma_plot_config
                    )

        self.current_iteration += 1
        if export_cplots:
            with suppress(TraitError):
                display(
                    widgets.VBox(
                        [
                            self.plasma_plot.show(
                                renderer="notebook_connected"
                            ),
                            self.luminosity_plot.show(
                                renderer="notebook_connected"
                            ),
                        ]
                    )
                )
