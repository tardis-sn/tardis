from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
from IPython.display import display
import matplotlib as mpl
import ipywidgets as widgets


def transistion_colors(name="jet", iterations=20):
    """
    Function to create colorscale for convergence plots, returns a list of colors.

    Parameters
    ----------
    name : string
        Name of the colorscale. Defaults to "jet".
    iterations : int
        Number of iterations. Defaults to 20.
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

        if "colorscale" in kwargs:
            self.colorscale = transistion_colors(
                kwargs["colorscale"], iterations=self.iterations
            )
        else:
            self.colorscale = transistion_colors()

        if "plasma_plot_config" in kwargs:
            if kwargs["plasma_plot_config"] != {}:
                self.plasma_plot_config = kwargs["plasma_plot_config"]

        if "luminosity_plot_config" in kwargs:
            if kwargs["luminosity_plot_config"] != {}:
                self.luminosity_plot_config = kwargs["luminosity_plot_config"]

    def fetch_data(self, name=None, value=None, type=None):
        """
        This allows user to fetch data from the Simulation class.
        This data is stored and used when an iteration is completed.

        Parameters
        ----------
        name : string
            name of the data
        value : string or array
            string or array of quantities,
        type : string
            either iterable or value

        """
        if type == "iterable":
            self.iterable_data[name] = value
        if type == "value":
            self.value_data[name].append(value)

    def create_plasma_plot(self):
        """
        Creates an empty plasma plot.
        The default layout can be overidden by passing plasma_plot_config dictionary in the run_tardis function.
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
            yaxis={"tickformat": "g", "title": r"$T_{rad}\ [K]$"},
            yaxis2={
                "tickformat": "g",
                "title": r"$W$",
                "range": [9000, 14000],
            },
            height=580,
        )

        # allows overriding default layout
        if hasattr(self, "plasma_plot_config"):
            for key in self.plasma_plot_config:
                fig["layout"][key] = self.plasma_plot_config[key]

        self.plasma_plot = fig

    def create_luminosity_plot(self):
        """
        Creates an empty luminosity plot.
        The default layout can be overidden by passing luminosity_plot_config dictionary in the run_tardis function.
        """
        marker_colors = ["#958aff", "#ff8b85", "#5cff74"]
        marker_line_colors = ["#27006b", "#800000", "#00801c"]
        marker_colors = ["#636EFA", "#EF553B", "#00CC96"]

        fig = go.FigureWidget().set_subplots(
            3,
            1,
            shared_xaxes=True,
            vertical_spacing=0.1,
            row_heights=[0.15, 0.7, 0.15],
        )

        for luminosity, marker_color, marker_line_color in zip(
            self.luminosities, marker_colors, marker_line_colors
        ):
            fig.add_scatter(
                name=luminosity + "<br>Luminosity",
                mode="lines+markers",
                row=2,
                col=1,
                marker_color=marker_color,
                marker_line_color=marker_line_color,
                legendgroup=luminosity,
                marker_line_width=1.5,
                opacity=0.6,
            )

        fig.add_scatter(
            name="Residual<br>Luminosity",
            row=3,
            col=1,
            hovertext="text",
            marker_color="rgb(158,202,225)",
            marker_line_color="rgb(8,48,107)",
            marker_line_width=1.5,
            mode="lines+markers",
            opacity=0.7,
        )

        fig.add_scatter(
            name="Inner<br>Boundary Temperature",
            row=1,
            col=1,
            hovertext="text",
            marker_color="rgb(158,202,225)",
            marker_line_color="rgb(8,48,107)",
            marker_line_width=1.5,
            mode="lines+markers",
            opacity=0.7,
        )

        fig = fig.update_layout(
            xaxis=dict(range=[0, self.iterations + 1], dtick=2),
            xaxis2=dict(
                matches="x",
                range=[0, self.iterations + 1],
                dtick=2,
            ),
            xaxis3=dict(
                matches="x",
                title=r"$\mbox{Iteration Number}$",
                dtick=2,
            ),
            yaxis=dict(
                title=r"$\mbox{T}_{inner}$",
                automargin=True,
                side="top",
                exponentformat="e",
            ),
            yaxis2=dict(
                exponentformat="e",
                title=r"$\mbox{Luminosity}~(erg~sec^{-1})$",
                title_font_size=13,
                automargin=True,
            ),
            yaxis3=dict(
                exponentformat="e",
                title=r"$~~\mbox{Residual}\\\mbox{Luminosity(%)}$",
                title_font_size=12,
                automargin=True,
            ),
            legend_tracegroupgap=0,
            hidesources=True,
            height=600,
            hoverlabel_align="right",
            legend_title_text="Luminosity",
        )

        # allows overriding default layout
        if hasattr(self, "luminosity_plot_config"):
            for key in self.luminosity_plot_config:
                fig["layout"][key] = self.luminosity_plot_config[key]

        self.luminosity_plot = fig

    def build(self):
        """
        Calls the create_plasma_plot and the create_luminosity_plot to build plots.
        """
        self.create_plasma_plot()
        self.create_luminosity_plot()
        display(
            widgets.VBox(
                [self.plasma_plot, self.luminosity_plot],
                layout=widgets.Layout(height="1000px"),
            )
        )

    def update_plasma_plots(self):
        """
        Updates the plasma plots using the data collected using the fetch_data function.
        This function is run every iteration.
        """
        x = self.iterable_data["velocity"].value.tolist()
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

        self.plasma_plot.add_scatter(
            x=[item / 100000 for item in x],
            y=self.iterable_data["t_rad"].value.tolist(),
            line_color=self.colorscale[self.current_iteration - 1],
            row=1,
            col=1,
            name=self.current_iteration,
            legendgroup=f"group-{self.current_iteration}",
            showlegend=False,
            customdata=customdata,
            hovertemplate="%{customdata}",
        )
        self.plasma_plot.add_scatter(
            x=[item / 100000 for item in x],
            y=self.iterable_data["w"].tolist(),  # TODO: is tolist() required?
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
        Updates the plasma plots using the data collected using the fetch_data function.
        This function is run every iteration.
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

    def update(self):
        """
        Calls functions used to build and update convergence plots.
        """
        if self.current_iteration == 1:
            self.build()
        self.update_plasma_plots()
        self.update_luminosity_plot()
        self.current_iteration += 1
