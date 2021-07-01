from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
from IPython.display import display
import matplotlib as mpl
import ipywidgets as widgets
from contextlib import suppress
from traitlets import TraitError
from astropy import units as u


def transition_colors(length, name="jet"):
    """
    Function to create colorscale for convergence plots, returns a list of colors.

    Parameters
    ----------
    length : int
        The length of the colorscale.
    name : string
        Name of the colorscale. Defaults to "jet".

    Returns
    -------
    colors: list
    """
    cmap = mpl.cm.get_cmap(name, length)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(mpl.colors.rgb2hex(rgb))
    return colors


class ConvergencePlots(object):
    """
    Class to create and update convergence plots for visualizing convergence of the \
    simulation.

    Parameters
    ----------
    iteration : int
        iteration number
    **kwargs : dict, optional
        keyword arguments

    Other Parameters
    ----------------
    plasma_plot_config, luminosity_plot_config : dict
        Dictionary used to override default plot properties of the plots.
        The  plasma_plot_config dictionary updates the plasma plot while 
        the luminosity_plot_config dict updates the luminosity and the inner boundary temperature plots.
        All properties related to data will be applied equally across all traces. 
    cmap : str, default: 'jet'
        String defining the cmap used in plots.
    export_cplots : bool, default: False
        If True, plots are displayed again using the `notebook_connected` renderer. This helps 
        display the plots in the documentation or in platforms like nbviewer. 
    
    """

    def __init__(self, iterations, **kwargs):
        self.iterable_data = {}
        self.value_data = defaultdict(list)
        self.iterations = iterations
        self.current_iteration = 1
        self.luminosities = ["Emitted", "Absorbed", "Requested"]
        self.plasma_plot = None
        self.luminosity_plot = None

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

        if "cmap" in kwargs:
            self.luminosity_line_colors = transition_colors(
                length=5,
                name=kwargs["cmap"],
            )
            self.plasma_colorscale = transition_colors(
                name=kwargs["cmap"], length=self.iterations
            )
        else:
            # default color scale is jet
            self.luminosity_line_colors = transition_colors(length=5)
            self.plasma_colorscale = transition_colors(length=self.iterations)

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
        # trace data for plasma plots is added in iterable data dictionary
        if item_type == "iterable":
            self.iterable_data[name] = value

        # trace data for luminosity plots and inner boundary temperature plot is stored in value_data dictionary
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
                "title": r"$\text{Velocity}~[\text{km}~\text{s}^{-1}]$",
            },
            xaxis2={
                "tickformat": "g",
                "title": r"$\text{Velocity}~[\text{km}~\text{s}^{-1}]$",
                "matches": "x",
            },
            yaxis={
                "tickformat": "g",
                "title": r"$T_{\text{rad}}\ [\text{K}]$",
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

        fig = go.FigureWidget().set_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            row_heights=[0.25, 0.5, 0.25],
        )

        fig.add_scatter(
            name="Inner<br>Boundary<br>Temperature",
            row=1,
            col=1,
            hovertext="text",
            marker_color=self.luminosity_line_colors[0],
            mode="lines",
        )

        for luminosity, line_color in zip(
            self.luminosities, self.luminosity_line_colors[1:4]
        ):
            fig.add_scatter(
                name=luminosity + "<br>Luminosity",
                mode="lines",
                row=2,
                col=1,
                marker_color=line_color,
            )

        fig.add_scatter(
            name="Residual<br>Luminosity",
            row=3,
            col=1,
            marker_color=self.luminosity_line_colors[4],
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
                title=r"$\mbox{Iteration Number}$",
                dtick=2,
            ),
            yaxis=dict(
                title=r"$T_{\text{inner}}\ [\text{K}]$",
                automargin=True,
                side="top",
                tickformat="g",
                exponentformat="e",
                tickmode="auto",
                nticks=4,
            ),
            yaxis2=dict(
                exponentformat="e",
                title=r"$\text{Luminosity}~[\text{erg s}^{-1}]$",
                title_font_size=13,
                automargin=True,
                tickmode="auto",
                nticks=7,
            ),
            yaxis3=dict(
                title=r"$~~\text{Residual}\\\text{Luminosity[%]}$",
                title_font_size=12,
                automargin=True,
                tickmode="auto",
                nticks=4,
            ),
            legend_tracegroupgap=0,
            height=630,
            hoverlabel_align="right",
            margin=dict(
                b=25, t=25, pad=0
            ),  # reduces whitespace surrounding the plot
        )

        # allows overriding default layout
        if hasattr(self, "luminosity_plot_config"):
            self.override_plot_parameters(fig, self.luminosity_plot_config)

        self.luminosity_plot = fig

    def override_plot_parameters(self, fig, parameters):
        """
        Overrides default plot properties. Any property inside the data dictionary is
        however, shared across all traces. This means trace-specific data properties can't be changed.

        Parameters
        ----------
        fig : go.FigureWidget
            FigureWidget object to be updated
        parameters : dict
            Dictionary used to update the default plot style. The dictionary should have a structure
            like that of go.FigureWidget.to_dict(), for more information please see https://plotly.com/python/figure-structure/
        """

        # because fig.data is a tuple of traces, a property in the data dictionary is applied to all traces
        # the fig is a nested dictionary, a property n levels deep is not changed until the value is a not dictionary
        # fig["property_1"]["property_2"]...["property_n"] = "value"
        for key, value in parameters.items():
            if key == "data":
                # all traces will have same data property
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
        # convert velocity to km/s
        x = self.iterable_data["velocity"].to(u.km / u.s).value.tolist()

        # add luminosity data in hover data in plasma plots
        customdata = len(x) * [
            "<br>"
            + "Emitted Luminosity: "
            + f'{self.value_data["Absorbed"][-1]:.4g}'
            + "<br>"
            + "Requested Luminosity: "
            + f'{self.value_data["Requested"][-1]:.4g}'
            + "<br>"
            + "Absorbed Luminosity: "
            + f'{self.value_data["Requested"][-1]:.4g}'
        ]

        # this adds radiation temperature vs shell velocity plot
        self.plasma_plot.add_scatter(
            x=x,
            y=self.iterable_data["t_rad"],
            line_color=self.plasma_colorscale[self.current_iteration - 1],
            row=1,
            col=1,
            name=self.current_iteration,
            legendgroup=f"group-{self.current_iteration}",
            showlegend=False,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.3f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

        # this adds dilution factor vs shell velocity plot
        self.plasma_plot.add_scatter(
            x=x,
            y=self.iterable_data["w"],
            line_color=self.plasma_colorscale[self.current_iteration - 1],
            row=1,
            col=2,
            legendgroup=f"group-{self.current_iteration}",
            name=self.current_iteration,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.3f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

    def update_luminosity_plot(self):
        """
        Updates luminosity convergence plots every iteration.
        """
        x = list(range(1, self.iterations + 1))

        with self.luminosity_plot.batch_update():
            # traces are updated according to the order they were added
            # the first trace is of the inner boundary temperature plot
            self.luminosity_plot.data[0].x = x
            self.luminosity_plot.data[0].y = self.value_data["t_inner"]
            self.luminosity_plot.data[
                0
            ].hovertemplate = "<b>%{y:.3f}</b> at X = %{x:,.0f}<extra>Inner Boundary Temperature</extra>"  # trace name in extra tag to avoid new lines in hoverdata

            # the next three for emitted, absorbed and requested luminosities
            for index, luminosity in zip(range(1, 4), self.luminosities):
                self.luminosity_plot.data[index].x = x
                self.luminosity_plot.data[index].y = self.value_data[luminosity]
                self.luminosity_plot.data[index].hovertemplate = (
                    "<b>%{y:.4g}</b>" + "<br>at X = %{x}<br>"
                )

                # last is for the residual luminosity
                y = [
                    ((emitted - requested) * 100) / requested
                    for emitted, requested in zip(
                        self.value_data["Emitted"], self.value_data["Requested"]
                    )
                ]

                self.luminosity_plot.data[4].x = x
                self.luminosity_plot.data[4].y = y
                self.luminosity_plot.data[
                    4
                ].hovertemplate = "<b>%{y:.2f}%</b> at X = %{x:,.0f}"

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
            # build only at first iteration
            if self.current_iteration == 1:
                self.build()

            self.update_plasma_plots()
            self.update_luminosity_plot()

            # data property for plasma plots needs to be
            # updated after the last iteration because new traces have been added
            if hasattr(self, "plasma_plot_config") and last:
                if "data" in self.plasma_plot_config:
                    self.override_plot_parameters(
                        self.plasma_plot, self.plasma_plot_config
                    )

        self.current_iteration += 1

        # the display function expects a Widget, while
        # fig.show() returns None, which causes the TraitError.
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
