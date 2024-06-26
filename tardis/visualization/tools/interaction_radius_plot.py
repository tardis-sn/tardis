import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import astropy.units as u

from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    species_string_to_tuple,
    roman_to_int,
    int_to_roman,
)
import tardis.visualization.tools.sdec_plot as sdec
from tardis.visualization import plot_util as pu


class InteractionRadiusPlotter:
    """
    Plotting interface for the interaction radius plot.
    """

    def __init__(self, data, time_explosion, no_of_shells):
        """
        Initialize the plotter with required data from the simulation.

        Parameters
        ----------
        data : dict of SDECData
            Dictionary to store data required for interaction radius plot,
            for both packet modes (real, virtual).

        time_explosion : astropy.units.Quantity
            Time of the explosion.

        no_of_shells : int
            The number of shells in TARDIS simulation.
        """

        self.data = data
        self.time_explosion = time_explosion
        self.velocity = velocity
        self.sdec_plotter = sdec.SDECPlotter(data)
        return

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of the plotter from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS simulation object produced by running a simulation.

        Returns
        -------
        InteractionRadiusPlotter
        """

        return cls(
            dict(
                virtual=sdec.SDECData.from_simulation(sim, "virtual"),
                real=sdec.SDECData.from_simulation(sim, "real"),
            ),
            sim.plasma.time_explosion,
            sim.simulation_state.no_of_shells,
        )

    @classmethod
    def from_hdf(cls, hdf_fpath):
        """
        Create an instance of the Plotter from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved.

        Returns
        -------
        InteractionRadiusPlotter
        """
        hdfstore = pd.HDFStore(hdf_fpath)
        time_explosion = (
            hdfstore["/simulation/plasma/scalars"]["time_explosion"] * u.s
        )
        velocity = (
            hdfstore["/simulation/simulation_state"]["velocity"] * u.cm / u.s
        )
        return cls(
            dict(
                virtual=sdec.SDECData.from_hdf(hdf_fpath, "virtual"),
                real=sdec.SDECData.from_hdf(hdf_fpath, "real"),
            ),
            time_explosion,
            velocity,
        )

    def _parse_species_list(self, species_list):
        """
        Parse user requested species list and create list of species ids to be used.

        Parameters
        ----------
        species_list : list of species to plot
            List of species (e.g. Si II, Ca II, etc.) that the user wants to show as unique colours.
            Species can be given as an ion (e.g. Si II), an element (e.g. Si), a range of ions
            (e.g. Si I - V), or any combination of these (e.g. species_list = [Si II, Fe I-V, Ca])

        Raises
        ------
        ValueError
            If species list contains invalid entries.

        """
        # Use the sdec.SDECPlotter instance to parse the species list
        self.sdec_plotter._parse_species_list(species_list)
        self._species_list = self.sdec_plotter._species_list
        self._species_mapped = self.sdec_plotter._species_mapped
        self._keep_colour = self.sdec_plotter._keep_colour

    def _make_colorbar_labels(self):
        """
        Generate labels for the colorbar based on species.

        If a species list is provided, uses that to generate labels.
        Otherwise, generates labels from the species in the model.
        """
        self.sdec_plotter.species = self.species
        self.sdec_plotter._make_colorbar_labels()
        self._species_name = self.sdec_plotter._species_name

    def _make_colorbar_colors(self):
        """Get the colours for the species to be plotted."""
        # the colours depends on the species present in the model and what's requested
        # some species need to be shown in the same colour, so the exact colours have to be
        # worked out

        color_list = []
        species_keys = list(self._species_mapped.keys())
        num_species = len(species_keys)

        for species_counter, species_key in enumerate(species_keys):
            if any(
                species in self.species
                for species in self._species_mapped[species_key]
            ):
                color = self.cmap(species_counter / num_species)
                color_list.append(color)

        self._color_list = color_list

    def _generate_plot_data(self, packets_mode):
        """
        Generate plot data and colors for species in the model.

        Parameters
        ----------
        packets_mode : str
            Packet mode, either 'virtual' or 'real'.

        Returns
        -------
        plot_data : list
            List of velocity data for each species.

        plot_colors : list
            List of colors corresponding to each species.
        """
        groups = self.data[packets_mode].packets_df_line_interaction.groupby(
            by="last_line_interaction_species"
        )

        plot_colors = []
        plot_data = []
        species_counter = 0

        for specie_list in self._species_mapped.values():
            full_v_last = []
            for specie in specie_list:
                if specie in self.species:
                    g_df = groups.get_group(specie)
                    r_last_interaction = (
                        g_df["last_interaction_in_r"].values * u.cm
                    )
                    v_last_interaction = (
                        r_last_interaction / self.time_explosion
                    ).to("km/s")
                    full_v_last.extend(v_last_interaction)
            if full_v_last:
                plot_data.append(full_v_last)
                plot_colors.append(self._color_list[species_counter])
                species_counter += 1

        return plot_data, plot_colors

    def _show_colorbar_mpl(self):
        """Show matplotlib colorbar with labels of elements mapped to colors."""

        if not self._species_name:
            raise ValueError("No species found to plot in colorbar")

        color_values = [
            self.cmap(species_counter / len(self._species_name))
            for species_counter in range(len(self._species_name))
        ]

        custcmap = clr.ListedColormap(color_values)
        norm = clr.Normalize(vmin=0, vmax=len(self._species_name))
        mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)
        mappable.set_array(np.linspace(1, len(self._species_name) + 1, 256))
        cbar = plt.colorbar(mappable, ax=self.ax)

        bounds = np.arange(len(self._species_name)) + 0.5
        cbar.set_ticks(bounds)

        cbar.set_ticklabels(self._species_name)
        return

    def generate_plot_mpl(
        self,
        packets_mode="virtual",
        ax=None,
        figsize=(11, 5),
        cmapname="jet",
        species_list=None,
    ):
        """
        Generate the last interaction radius distribution plot using matplotlib.

        Parameters
        ----------
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, creates a new figure.
        figsize : tuple, optional
            Size of the figure. Default is (11, 6).
        cmapname : str, optional
            Colormap name. Default is 'jet'. A specific colormap can be chosen, such as "jet", "viridis", "plasma", etc.
        species_list : list of str
            List of species to plot.

        Returns
        -------
        matplotlib.axes.Axes
            Axes object with the plot.
        """

        # Parse the requested species list
        self._parse_species_list(species_list)
        species_in_model = np.unique(
            self.data[packets_mode]
            .packets_df_line_interaction["last_line_interaction_species"]
            .values
        )
        if self._species_list is None or not self._species_list:
            raise ValueError("No species provided for plotting.")
        msk = np.isin(self._species_list, species_in_model)
        self.species = np.array(self._species_list)[msk]

        if len(self.species) == 0:
            raise ValueError("No valid species found for plotting.")

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        # Get the labels in the color bar. This determines the number of unique colors
        self._make_colorbar_labels()
        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = cm.get_cmap(cmapname, len(self._species_name))
        # Get the number of unique colors
        self._make_colorbar_colors()
        # Show colorbar
        # self._show_colorbar_mpl()

        plot_data, plot_colors = self._generate_plot_data(packets_mode)
        bin_edges = (self.velocity).to("km/s")

        for data, color, name in zip(
            plot_data, plot_colors, self._species_name
        ):
            hist, _ = np.histogram(data, bins=bin_edges)
            step_x = np.repeat(bin_edges, 2)[1:-1]
            step_y = np.repeat(hist, 2)
            self.ax.plot(
                step_x,
                step_y,
                label=name,
                color=color,
                linewidth=2.5,
                drawstyle="steps-post",
                alpha=0.75,
            )

        self.ax.ticklabel_format(axis="y", scilimits=(0, 0))
        self.ax.tick_params("both", labelsize=14)
        self.ax.set_xlabel("Last Interaction Velocity (km/s)", fontsize=14)
        self.ax.set_ylabel("Packet Count", fontsize=14)
        self.ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        self.ax.legend(fontsize=15)
        plt.tight_layout()

        return self.ax

    # def _show_colorbar_mpl(self):
    #     """Show matplotlib colorbar with labels of elements mapped to colors."""

    #     color_values = [
    #         self.cmap(species_counter / len(self._species_name))
    #         for species_counter in range(len(self._species_name))
    #     ]

    #     custcmap = clr.ListedColormap(color_values)
    #     norm = clr.Normalize(vmin=0, vmax=len(self._species_name))
    #     mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)
    #     mappable.set_array(np.linspace(1, len(self._species_name) + 1, 256))
    #     cbar = plt.colorbar(mappable, ax=self.ax)

    #     bounds = np.arange(len(self._species_name)) + 0.5
    #     cbar.set_ticks(bounds)

    #     cbar.set_ticklabels(self._species_name)
    #     return

    def generate_plot_ply(
        self,
        packets_mode="virtual",
        species_list=None,
        cmapname="jet",
    ):
        """
        Generate the last interaction radius distribution plot using plotly.

        Parameters
        ----------
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        species_list : list of str
            List of species to plot.

        Returns
        -------
        plotly.graph_objects.Figure
            Plotly figure object with the plot.
        """
        # Parse the requested species list
        self._parse_species_list(species_list)
        species_in_model = np.unique(
            self.data[packets_mode]
            .packets_df_line_interaction["last_line_interaction_species"]
            .values
        )
        if self._species_list is None or not self._species_list:
            raise ValueError("No species provided for plotting.")
        msk = np.isin(self._species_list, species_in_model)
        self.species = np.array(self._species_list)[msk]

        if len(self.species) == 0:
            raise ValueError("No valid species found for plotting.")
        # Get the labels in the color bar. This determines the number of unique colors
        self._make_colorbar_labels()
        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = cm.get_cmap(cmapname, len(self._species_name))
        # Get the number of unique colors
        self._make_colorbar_colors()

        plot_data, plot_colors = self._generate_plot_data(packets_mode)
        bin_edges = (self.velocity).to("km/s")
        fig = go.Figure()
        for data, color, name in zip(
            plot_data, plot_colors, self._species_name
        ):
            hist, _ = np.histogram(data, bins=bin_edges)
            self.step_x = np.repeat(bin_edges, 2)[1:-1]
            self.step_y = np.repeat(hist, 2)
            fig.add_trace(
                go.Scatter(
                    x=self.step_x,
                    y=self.step_y,
                    mode="lines",
                    line=dict(
                        color=pu.to_rgb255_string(color),
                        width=2.5,
                        shape="hv",
                    ),
                    name=name,
                    opacity=0.75,
                    nbinsx=self.no_of_shells,
                )
            )
        fig.update_layout(
            barmode="overlay",
            xaxis_title="Last Interaction Velocity (km/s)",
            yaxis_title="Packet Count",
            font=dict(size=14),
            yaxis=dict(tickformat=".1e"),
            xaxis=dict(tickformat=".0f", tickmode="auto"),
        )

        # fig = self._show_colorbar_ply(fig)
        return fig

    # def _show_colorbar_ply(self, fig):
    #     """
    #     Show plotly colorbar with labels of elements mapped to colors.

    #     Parameters
    #     ----------
    #     fig : plotly.graph_objects.Figure
    #     Plotly figure object to add the colorbar to.

    #     Returns
    #     -------
    #     plotly.graph_objects.Figure
    #         Plotly figure object with the colorbar added.
    #     """
    #     # Interpolate [0, 1] range to create bins equal to number of elements
    #     colorscale_bins = np.linspace(0, 1, num=len(self._species_name) + 1)

    #     # Create a categorical colorscale [a list of (reference point, color)]
    #     # by mapping same reference points (excluding 1st and last bin edge)
    #     # twice in a row (https://plotly.com/python/colorscales/#constructing-a-discrete-or-discontinuous-color-scale)
    #     categorical_colorscale = []
    #     for species_counter in range(len(self._species_name)):
    #         color = pu.to_rgb255_string(
    #             self.cmap(colorscale_bins[species_counter])
    #         )
    #         categorical_colorscale.append(
    #             (colorscale_bins[species_counter], color)
    #         )
    #         categorical_colorscale.append(
    #             (colorscale_bins[species_counter + 1], color)
    #         )

    #     coloraxis_options = {
    #         "colorscale": categorical_colorscale,
    #         "showscale": True,
    #         "cmin": 0,
    #         "cmax": len(self._species_name),
    #         "colorbar": {
    #             "title": "Elements",
    #             "tickvals": np.arange(0, len(self._species_name)) + 0.5,
    #             "ticktext": self._species_name,
    #             # to change length and position of colorbar
    #             "len": 1,
    #             "yanchor": "top",
    #             "y": 1,
    #         },
    #     }

    #     colorbar_trace = go.Scatter(
    #         x=[None],
    #         y=[0],
    #         mode="markers",
    #         name="Colorbar",
    #         showlegend=False,
    #         hoverinfo="skip",
    #         marker=dict(color=[0], opacity=0, **coloraxis_options),
    #     )

    #     fig.add_trace(colorbar_trace)
    #     return fig
