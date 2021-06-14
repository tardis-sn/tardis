from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
import matplotlib as mpl
import ipywidgets as widgets


def transistion_colors(name="jet", iterations=20):
    cmap = mpl.cm.get_cmap(name, iterations)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(mpl.colors.rgb2hex(rgb))
    return colors


class ConvergencePlots:
    def __init__(self):
        self.iterable_data = {}
        self.value_data = defaultdict(list)
        self.rows = 4
        self.cols = 2
        self.specs = [
            [{}, {}],
            [{"colspan": 2}, None],
            [{"colspan": 2}, None],
            [{"colspan": 2}, None],
        ]
        self.row_heights = [0.45, 0.1, 0.4, 0.1]
        self.vertical_spacing = 0.07
        self.iterations = 20
        self.current_iteration = 1

    def fetch_data(self, name=None, value=None):
        """
        This allows user to fetch data from the Simulation class.
        This data is stored and used when an iteration is completed.
        Returns:
            iterable_data: iterable data from the simulation, values like radiation temperature
            value_data: values like luminosity
        """
        self.iterable_data[name] = value
        self.value_data[name].append(value)
        self.current_iteration += 1

    def get_data(self):
        """
        Returns data for convergence plots

        Returns:
            iterable_data: iterable data from the simulation, values like radiation temperature
            value_data: values like luminosity
        """
        return self.iterable_data, self.value_data


class BuildCplots(ConvergencePlots):
    def __init__(self):
        super().__init__()
        self.use_vbox = True

    def create_plasma_plot(self):
        """
        creates empty plasma plot
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
                "range": [9000, 14000],
            },
            yaxis2={"tickformat": "g", "title": r"$W$"},
            height=580,
        )

        self.plasma_plot = fig

    def create_luminosity_plot(self):
        marker_colors = ["#958aff", "#ff8b85", "#5cff74"]
        marker_line_colors = ["#27006b", "#800000", "#00801c"]
        marker_colors = ["#636EFA", "#EF553B", "#00CC96"]
        self.luminosities = ["Emitted", "Absorbed", "Requested"]

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
            name="Next Inner<br>Boundary Temperature",
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
            xaxis=dict(range=[0, 21], dtick=2),
            xaxis2=dict(
                matches="x",
                range=[0, 21],
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
        self.luminosity_plot = fig

    def build(self):
        self.create_plasma_plot()
        self.create_luminosity_plot()

        return widgets.VBox([self.plasma_plot, self.luminosity_plot])


class UpdateCplots(BuildCplots):
    def __init__(self):
        super().__init__()

    def update_plasma_plots(self):
        x = self.iterable_data["velocity"].value.tolist()
        # customdata = len(shell_velocity)*['Emitted Luminosity: ' +f'{emitted_luminosity.value:.2g}' +'<br>'+'Requested Luminosity: ' + f'{reabsorbed_luminosity.value:.2g}' + '<br>'+'Absorbed Luminosity: ' + f'{sim.luminosity_requested.value:.2g}']

        self.plasma_plot.add_scatter(
            x=[item / 100000 for item in x],
            y=self.iterable_data["t_rad"].value.tolist(),
            line_color=transistion_colors()[self.current_iteration - 1],
            row=1,
            col=1,
            name=self.current_iteration,
            legendgroup=f"group-{self.current_iteration}",
            showlegend=False,
            # customdata = customdata,
            hovertemplate="%{customdata}",
        )
        self.plasma_plot.add_scatter(
            x=[item / 100000 for item in x],
            y=self.iterable_data[
                "w"
            ].value.tolist(),  # TODO: is tolist() required?
            line_color=transistion_colors()[self.current_iteration - 1],
            row=1,
            col=2,
            legendgroup=f"group-{self.current_iteration}",
            name=self.current_iteration,
            # customdata = customdata,
            hovertemplate="%{customdata}",
        )

    def update_luminosity_plot(self):
        x = list(range(1, self.iterations + 1))
        with self.luminosity_plot.batch_update():
            for index, luminosity in zip(range(3), self.luminosities):
                self.luminosity_plot.data[index].x = x
                self.luminosity_plot.data[index].y = self.value_data[
                    luminosity
                ].value
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
                ].hovertemplate = "Next Inner Body Temperature: %{y:.2f} at X = %{x:,.0f}<extra></extra>"
