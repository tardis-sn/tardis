import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
from pylab import *  # TODO: fix this

i, ic = 0, 0


def transistion_colors(name="Blues", iterations=20):
    cmap = cm.get_cmap(name, iterations)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(matplotlib.colors.rgb2hex(rgb))
    return colors


class build_cplots:
    def __init__(self):
        # TODO: ask for parameters
        self.plasma_shell_v_plot = None

    # TODO:  more convergence plots functions
    def plasma_shell_v(self):
        fig = go.FigureWidget().set_subplots(1, 2, shared_xaxes=True)
        fig.add_scatter(row=1, col=1, legendgroup="group")
        fig.add_scatter(row=1, col=2, legendgroup="group")

        fig = fig.update_layout(
            showlegend=True,
            yaxis=dict(
                title=r"$W$",
            ),
            yaxis2=dict(title=r"$T_{rad}\ [K]$", range=[9000, 14000]),
            xaxis=dict(
                showexponent="all",
                title=r"$Shell~~Velocity$",
                matches="x2",
                exponentformat="e",
            ),
            xaxis2=dict(
                visible=True,
                showexponent="all",
                title=r"$Shell~~Velocity$",
                exponentformat="e",
            ),
        )
        self.plasma_shell_v_plot = fig
        return self.plasma_shell_v_plot


class update_cplots:
    def __init__(self, w, t_rad, velocity, plasma_shell_v_plot=None):
        self.w = w
        self.t_rad = t_rad  # list
        self.velocity = velocity  # list
        self.plasma_shell_v_plot = plasma_shell_v_plot
        self.colors = transistion_colors()

    def update_plasma_shell_v_colors(self, index):
        self.plasma_shell_v_plot["data"][index]["line"]["color"] = self.colors[
            ic - 1
        ]  # TODO: ic

    def update_plasma_shell_v_data(self):
        # update t_rad subplot
        self.plasma_shell_v_plot.add_scatter(
            x=self.velocity,
            y=self.t_rad,
            line_color=self.colors[ic],
            row=1,
            col=2,
            # name = f'iteration-{ic}', # TODO
            # legendgroup=f"group-{ic}",
            showlegend=False,
        )

        self.plasma_shell_v_plot.add_scatter(
            x=self.velocity,
            y=self.w,
            line_color=self.colors[ic],
            row=1,
            col=2,
            # name = f'iteration-{ic}', # TODO
            # legendgroup=f"group-{ic}",
            showlegend=False,
        )
