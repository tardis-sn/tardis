"""Convergence Plots to see the convergence of the simulation in real time."""

from collections import defaultdict
import warnings
import os
import csv

import ipywidgets as widgets
import numpy as np
from astropy import units as u
from IPython.display import HTML, display

# Added the below as a (temporary) workaround to the latex
# labels on the convergence plots not rendering correctly.
import plotly
import plotly.graph_objects as go
import plotly.io as pio

import tardis.visualization.plot_util as pu

plotly.offline.init_notebook_mode(connected=True)
# mathjax needs to be loaded for latex labels to render correctly
# see https://github.com/tardis-sn/tardis/issues/2446
display(
    HTML(
        '<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_SVG"></script>'
    )
)


class ConvergencePlots:
    """
    Create and update convergence plots for visualizing convergence of the simulation.

    Parameters
    ----------
    iterations : int
        iteration number
    **kwargs : dict, optional
        Additional keyword arguments. These arguments are defined in the Other Parameters section.

    Other Parameters
    ----------------
    plasma_plot_config : dict, optional
        Dictionary used to override default plot properties of plasma plots.
    t_inner_luminosities_config : dict, optional
        Dictionary used to override default plot properties of the inner boundary temperature and luminosity plots.
    plasma_cmap : str, default: 'jet', optional
        String defining the cmap used in plasma plots.
    t_inner_luminosities_colors : str or list, optional
        String defining cmap for luminosity and inner boundary temperature plot.
        The list can be a list of colors in rgb, hex or css-names format as well.
    export_convergence_plots : bool, default: False, optional
        If True, plots are displayed again using the `notebook_connected` renderer. This helps
        to display the plots in the documentation or in platforms like nbviewer.
    export_format : str, optional
        If provided, file export is performed on the last iteration in the specified format.
        Supported formats: 'ascii', 'csv' for data export, or image formats supported by
        plotly/kaleido (e.g. 'png', 'svg', 'pdf'). File export is independent of export_convergence_plots.
    export_path : str, optional
        Directory where exported files will be written. Defaults to the current working
        directory.
    export_prefix : str, optional
        Filename prefix used for exported files. Defaults to 'convergence_plots'.

    Notes
    -----
        When overriding plot's configuration using the `plasma_plot_config` and the
        `t_inner_luminosities_config` dictionaries, data related properties are
        applied equally accross all traces.
        The dictionary should have a structure like that of `plotly.graph_objs.FigureWidget.to_dict()`,
        for more information please see https://plotly.com/python/figure-structure/
    """

    def __init__(self, iterations, **kwargs):
        self.iterable_data = {}
        self.value_data = defaultdict(list)
        self.iterations = iterations
        self.current_iteration = 1
        self.luminosities = ["Emitted", "Absorbed", "Requested"]
        self.plasma_plot = None
        self.t_inner_luminosities_plot = None
        self.display_handle = None
        self.export_format = kwargs.get("export_format", None)
        self.export_path = kwargs.get("export_path", None)
        self.export_prefix = kwargs.get("export_prefix", "convergence_plots")

        if "plasma_plot_config" in kwargs:
            self.plasma_plot_config = kwargs["plasma_plot_config"]

        if "t_inner_luminosities_config" in kwargs:
            self.t_inner_luminosities_config = kwargs["t_inner_luminosities_config"]

        if "plasma_cmap" in kwargs:
            self.plasma_colorscale = pu.get_hex_color_strings(
                name=kwargs["plasma_cmap"], length=self.iterations
            )
        else:
            # default color scale is jet
            self.plasma_colorscale = pu.get_hex_color_strings(length=self.iterations)

        if "t_inner_luminosities_colors" in kwargs:
            # use cmap if string
            if type(kwargs["t_inner_luminosities_colors"]) == str:
                self.t_inner_luminosities_colors = pu.get_hex_color_strings(
                    length=5,
                    name=kwargs["t_inner_luminosities_colors"],
                )
            else:
                self.t_inner_luminosities_colors = kwargs["t_inner_luminosities_colors"]
        else:
            # using default plotly colors
            self.t_inner_luminosities_colors = [None] * 5

    def fetch_data(self, name=None, value=None, item_type=None):
        """
        Fetch data from the Simulation class.

        Parameters
        ----------
        name : string
            name of the data
        value : string or array
            string or an array of quantities
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
        """Create an empty plasma plot."""
        fig = go.FigureWidget().set_subplots(rows=1, cols=2, shared_xaxes=True)

        # empty traces to build figure
        fig.add_scatter(row=1, col=1)
        fig.add_scatter(row=1, col=2)

        # 2 y axes and 2 x axes correspond to the 2 subplots in the plasma plot
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
                "nticks": 15,
            },
            yaxis2={
                "tickformat": "g",
                "title": r"$W$",
                "nticks": 15,
            },
            height=450,
            legend_title_text="Iterations",
            legend_traceorder="reversed",
            margin=dict(
                l=10, r=135, b=25, t=25, pad=0
            ),  # reduce whitespace surrounding the plot and increase right indentation to align with the t_inner and luminosity plot
        )

        # allow overriding default layout
        if hasattr(self, "plasma_plot_config"):
            self.override_plot_parameters(fig, self.plasma_plot_config)
        self.plasma_plot = fig

    def create_t_inner_luminosities_plot(self):
        """Create an empty t_inner and luminosity plot."""
        fig = go.FigureWidget().set_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            row_heights=[0.25, 0.5, 0.25],
        )

        # add inner boundary temperature vs iterations plot
        fig.add_scatter(
            name="Inner<br>Boundary<br>Temperature",
            row=1,
            col=1,
            hovertext="text",
            marker_color=self.t_inner_luminosities_colors[0],
            mode="lines",
        )

        # add luminosity vs iterations plot
        # has three traces for emitted, requested and absorbed luminosities
        for luminosity, line_color in zip(
            self.luminosities, self.t_inner_luminosities_colors[1:4]
        ):
            fig.add_scatter(
                name=luminosity + "<br>Luminosity",
                mode="lines",
                row=2,
                col=1,
                marker_color=line_color,
            )

        # add residual luminosity vs iterations plot
        fig.add_scatter(
            name="Residual<br>Luminosity",
            row=3,
            col=1,
            marker_color=self.t_inner_luminosities_colors[4],
            mode="lines",
        )

        # 3 y axes and 3 x axes correspond to the 3 subplots in the t_inner and luminosity convergence plot
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
                tickformat="g",
                exponentformat="e",
                nticks=4,
            ),
            yaxis2=dict(
                exponentformat="e",
                title=r"$\text{Luminosity}~[\text{erg s}^{-1}]$",
                title_font_size=13,
                automargin=True,
                nticks=7,
            ),
            yaxis3=dict(
                title=r"$~~\text{Residual}\\\text{Luminosity[%]}$",
                title_font_size=12,
                automargin=True,
                nticks=4,
            ),
            height=630,
            hoverlabel_align="right",
            margin=dict(b=25, t=25, pad=0),  # reduces whitespace surrounding the plot
        )

        # allow overriding default layout
        if hasattr(self, "t_inner_luminosities_config"):
            self.override_plot_parameters(fig, self.t_inner_luminosities_config)

        self.t_inner_luminosities_plot = fig

    def override_plot_parameters(self, fig, parameters):
        """
        Override default plot properties.

        Any property inside the data dictionary is however, applied equally across all traces.
        This means trace-specific data properties can't be changed using this function.

        Parameters
        ----------
        fig : go.FigureWidget
            FigureWidget object to be updated
        parameters : dict
            Dictionary used to update the default plot style.
        """
        # because fig.data is a tuple of traces, a property in the data dictionary is applied to all traces
        # the fig is a nested dictionary, any property n levels deep is not changed until the value is a not dictionary
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
        Create empty convergence plots and display them.

        Parameters
        ----------
        display_plot : bool, default: True, optional
            Displays empty plots.
        """
        self.create_plasma_plot()
        self.create_t_inner_luminosities_plot()

        if display_plot:
            self.display_handle = display(
                widgets.VBox(
                    [self.plasma_plot, self.t_inner_luminosities_plot],
                ),
                display_id="convergence_plots",
            )

    def update_plasma_plots(self):
        """Update plasma convergence plots every iteration."""
        # convert velocity to km/s
        velocity_km_s = self.iterable_data["velocity"].to(u.km / u.s).value.tolist()

        # add luminosity data in hover data in plasma plots
        customdata = len(velocity_km_s) * [
            "<br>"
            "Emitted Luminosity: "
            f"{self.value_data['Emitted'][-1]:.4g}"
            "<br>"
            "Requested Luminosity: "
            f"{self.value_data['Requested'][-1]:.4g}"
            "<br>"
            "Absorbed Luminosity: "
            f"{self.value_data['Absorbed'][-1]:.4g}"
        ]

        # add a radiation temperature vs shell velocity trace to the plasma plot
        self.plasma_plot.add_scatter(
            x=velocity_km_s,
            y=np.append(self.iterable_data["t_rad"], self.iterable_data["t_rad"][-1:]),
            line_color=self.plasma_colorscale[self.current_iteration - 1],
            line_shape="hv",
            row=1,
            col=1,
            name=self.current_iteration,
            legendgroup=f"group-{self.current_iteration}",
            showlegend=False,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.3f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

        # add a dilution factor vs shell velocity trace to the plasma plot
        self.plasma_plot.add_scatter(
            x=velocity_km_s,
            y=np.append(self.iterable_data["w"], self.iterable_data["w"][-1:]),
            line_color=self.plasma_colorscale[self.current_iteration - 1],
            line_shape="hv",
            row=1,
            col=2,
            legendgroup=f"group-{self.current_iteration}",
            name=self.current_iteration,
            customdata=customdata,
            hovertemplate="<b>Y</b>: %{y:.3f} at <b>X</b> = %{x:,.0f}%{customdata}",
        )

    def update_t_inner_luminosities_plot(self):
        """Update the t_inner and luminosity convergence plots every iteration."""
        x = list(range(1, self.iterations + 1))

        with self.t_inner_luminosities_plot.batch_update():
            # traces are updated according to the order they were added
            # the first trace is of the inner boundary temperature plot
            self.t_inner_luminosities_plot.data[0].x = x
            self.t_inner_luminosities_plot.data[0].y = self.value_data["t_inner"]
            self.t_inner_luminosities_plot.data[0].hovertemplate = (
                "<b>%{y:.3f}</b> at X = %{x:,.0f}<extra>Inner Boundary Temperature</extra>"  # trace name in extra tag to avoid new lines in hoverdata
            )

            # the next three for emitted, absorbed and requested luminosities
            for index, luminosity in zip(range(1, 4), self.luminosities):
                self.t_inner_luminosities_plot.data[index].x = x
                self.t_inner_luminosities_plot.data[index].y = self.value_data[
                    luminosity
                ]
                self.t_inner_luminosities_plot.data[index].hovertemplate = (
                    "<b>%{y:.4g}</b>" + "<br>at X = %{x}<br>"
                )

            # last is for the residual luminosity
            y = [
                ((emitted - requested) * 100) / requested
                for emitted, requested in zip(
                    self.value_data["Emitted"], self.value_data["Requested"]
                )
            ]

            self.t_inner_luminosities_plot.data[4].x = x
            self.t_inner_luminosities_plot.data[4].y = y
            self.t_inner_luminosities_plot.data[4].hovertemplate = (
                "<b>%{y:.2f}%</b> at X = %{x:,.0f}"
            )

    def update(self, export_convergence_plots=False, last=False, export_format=None, export_path=None, export_prefix="convergence_plots"):
        """
        Update the convergence plots every iteration.

        Parameters
        ----------
        export_convergence_plots : bool, default: False, optional
            Displays the convergence plots again using plotly's `notebook_connected` renderer.
            This helps to display the plots in notebooks when shared on platforms like nbviewer.
            Please see https://plotly.com/python/renderers/ for more information.
        last : bool, default: False, optional
            True if it's last iteration.
        export_format : str, optional
            If provided, file export is triggered on the last iteration in the specified format.
            Supported formats: 'ascii', 'csv' for data export, or image formats supported by
            plotly/kaleido (e.g. 'png', 'svg', 'pdf'). File export is independent of export_convergence_plots.
        export_path : str, optional
            Directory where exported files will be written. Defaults to the current working
            directory.
        export_prefix : str, default: 'convergence_plots', optional
            Filename prefix used for exported files.
        """
        if self.iterable_data != {}:
            # build only at first iteration
            if self.current_iteration == 1:
                self.build()

            self.update_plasma_plots()
            self.update_t_inner_luminosities_plot()

            # Update the existing widget display unless we're in export mode
            # In export mode, we'll display the plots using show() instead
            if self.display_handle is not None and not (
                export_convergence_plots and last
            ):
                self.display_handle.update(
                    widgets.VBox(
                        [self.plasma_plot, self.t_inner_luminosities_plot],
                    )
                )

            # data property for plasma plots needs to be
            # updated after the last iteration because new traces have been added
            if hasattr(self, "plasma_plot_config") and last:
                if "data" in self.plasma_plot_config:
                    self.override_plot_parameters(
                        self.plasma_plot, self.plasma_plot_config
                    )

        self.current_iteration += 1

        # Resolve effective export parameters (allow instance-level defaults)
        efmt = export_format or getattr(self, "export_format", None)
        epath = export_path or getattr(self, "export_path", os.getcwd())
        prefix = export_prefix or getattr(self, "export_prefix", "convergence_plots")

        # Triggered on last=True when export_format is set
        if last and efmt and self.plasma_plot is not None and self.t_inner_luminosities_plot is not None:
            try:
                os.makedirs(epath, exist_ok=True)
            except OSError:
                warnings.warn(f"Could not create export directory '{epath}'")

            if efmt in ("ascii", "csv"):
                try:
                    self._export_csv(epath, prefix)
                except (OSError, csv.Error) as e:
                    warnings.warn(f"Failed to export CSV data to {epath}. Error: {e}")
            else:
                try:
                    self._export_images(epath, prefix, efmt)
                except (OSError, ValueError, ImportError) as e:
                    warnings.warn(f"Failed to export images to {epath}. Error: {e}")

        # Export convergence plots for sharing on platforms like nbviewer
        # This creates additional outputs that work better for static notebook sharing
        # but can cause duplicate displays in interactive environments
        if (
            export_convergence_plots
            and last
            and self.plasma_plot is not None
            and self.t_inner_luminosities_plot is not None
        ):
            self.plasma_plot.show(renderer="notebook_connected")
            self.t_inner_luminosities_plot.show(renderer="notebook_connected")

    def _export_images(self, export_path, prefix, fmt):
        """
        Export plots as static images using Plotly + Kaleido.

        Parameters
        ----------
        export_path : str
            Directory where images will be written.
        prefix : str
            Filename prefix.
        fmt : str
            Image format (e.g. 'png', 'svg', 'pdf').
        """
        plasma_file = os.path.join(export_path, f"{prefix}_plasma.{fmt}")
        tinner_file = os.path.join(export_path, f"{prefix}_tinner.{fmt}")

        try:
            plasma_fig = go.Figure(self.plasma_plot.to_plotly_json())
            tinner_fig = go.Figure(self.t_inner_luminosities_plot.to_plotly_json())
        except ValueError:
            plasma_fig = go.Figure(self.plasma_plot)
            tinner_fig = go.Figure(self.t_inner_luminosities_plot)

        try:
            import kaleido
        except ImportError:
            warnings.warn(
                "Kaleido is not available; image export may fail. Install 'kaleido' to enable image export."
            )

        try:
            plasma_fig.write_image(plasma_file)
            tinner_fig.write_image(tinner_file)
        except (ValueError, OSError, ImportError) as e:
            warnings.warn(
                f"Failed to write images to {export_path}. Ensure 'kaleido' is installed. Error: {e}"
            )

    def _export_csv(self, export_path, prefix):
        """
        Export convergence history and last-iteration plasma data as CSV/ASCII.

        Parameters
        ----------
        export_path : str
            Directory where CSV files will be written.
        prefix : str
            Filename prefix.
        """

        tinner_file = os.path.join(export_path, f"{prefix}_tinner.csv")
        tinner_values = self.value_data.get("t_inner", [])
        emitted_values = self.value_data.get("Emitted", [])
        absorbed_values = self.value_data.get("Absorbed", [])
        requested_values = self.value_data.get("Requested", [])

        iterations = list(range(1, len(tinner_values) + 1))

        residuals = []
        for e, r in zip(emitted_values, requested_values):
            if e is None or r is None:
                residuals.append(None)
                continue
            r_val = r.value if hasattr(r, "value") else r
            if r_val == 0:
                residuals.append(None)
            else:
                residuals.append(((e - r) * 100) / r)

        with open(tinner_file, "w", newline="") as fh:
            writer = csv.writer(fh)
            header = [
                "iteration",
                "t_inner [K]",
                "Emitted [erg s^-1]",
                "Absorbed [erg s^-1]",
                "Requested [erg s^-1]",
                "Residual_percent [%]",
            ]
            writer.writerow(header)

            for i in range(len(tinner_values)):
                tval = tinner_values[i]
                if hasattr(tval, "to") and tval.unit.is_equivalent(u.K):
                    tval_out = float(tval.to(u.K).value)
                else:
                    tval_out = float(getattr(tval, "value", tval))

                def _val_at(arr, idx):
                    if idx >= len(arr):
                        return None
                    v = arr[idx]
                    if hasattr(v, "to") and v.unit.is_equivalent(u.erg / u.s):
                        return float(v.to(u.erg / u.s).value)
                    return float(getattr(v, "value", v))

                row = [
                    iterations[i],
                    tval_out,
                    _val_at(emitted_values, i),
                    _val_at(absorbed_values, i),
                    _val_at(requested_values, i),
                    residuals[i] if i < len(residuals) else None,
                ]
                writer.writerow(row)

        plasma_file = os.path.join(export_path, f"{prefix}_plasma_last.csv")
        vel = self.iterable_data.get("velocity")
        t_rad = self.iterable_data.get("t_rad")
        w = self.iterable_data.get("w")

        if vel is None or t_rad is None or w is None:
            warnings.warn("Plasma data not available for CSV export; skipping plasma CSV.")
            return

        if hasattr(vel, "to") and vel.unit.is_equivalent(u.km / u.s):
            vel_vals = list(vel.to(u.km / u.s).value)
        else:
            vel_vals = [
                float(v.to(u.km / u.s).value) if hasattr(v, "to") and v.unit.is_equivalent(u.km / u.s) else float(getattr(v, "value", v))
                for v in vel
            ]

        if hasattr(t_rad, "to") and t_rad.unit.is_equivalent(u.K):
            tvals = list(t_rad.to(u.K).value)
        else:
            tvals = [
                float(t.to(u.K).value) if hasattr(t, "to") and t.unit.is_equivalent(u.K) else float(getattr(t, "value", t))
                for t in t_rad
            ]

        w_vals = [float(getattr(wi, "value", wi)) for wi in w]

        with open(plasma_file, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["velocity [km s^-1]", "t_rad [K]", "w [dimensionless]"])
            for xv, yv, wv in zip(vel_vals, tvals, w_vals):
                writer.writerow([xv, yv, wv])
