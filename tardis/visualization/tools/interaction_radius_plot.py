import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import astropy.units as u

from tardis.util.base import (
    atomic_number2element_symbol,
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
        if self._species_list is None:
            species_name = [
                atomic_number2element_symbol(atomic_num)
                for atomic_num in self.species
            ]
        else:
            species_name = []
            for species_key, species_ids in self._species_mapped.items():
                if any(species in self.species for species in species_ids):
                    if species_key % 100 == 0:
                        label = atomic_number2element_symbol(species_key // 100)
                    else:
                        atomic_number = species_key // 100
                        ion_number = species_key % 100
                        ion_numeral = int_to_roman(ion_number + 1)
                        label = f"{atomic_number2element_symbol(atomic_number)} {ion_numeral}"
                    species_name.append(label)

        self._species_name = species_name

    def _make_colorbar_colors(self):
        """
        Generate colors for the species to be plotted.

        Creates a list of colors corresponding to the species names.
        """
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

    def _prepare_plot(self, packets_mode, species_list, cmapname):
        """
        Helper function to prepare the plot data and colors.

        Parameters
        ----------
        packets_mode : str
            Packet mode, either 'virtual' or 'real'.
        species_list : list of str
            List of species to plot.

        Raises
        ------
        ValueError
            If no species provided or no valid species found.

        """
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

        self._make_colorbar_labels()
        self.cmap = cm.get_cmap(cmapname, len(self._species_name))
        self._make_colorbar_colors()

    def generate_plot_mpl(
        self,
        packets_mode="virtual",
        ax=None,
        figsize=(11, 5),
        cmapname="jet",
        species_list=None,
        log_scale=False,
        bins_range=None,
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
        log_scale : bool, optional
            If True, both axes are scaled logarithmically. Default is False.
        bins_range : tuple of int, optional
            Range of bins to plot, e.g., (3, 9). Default is None, which plots all bins.

        Returns
        -------
        matplotlib.axes.Axes
            Axes object with the plot.
        """
        self._prepare_plot(packets_mode, species_list, cmapname)
        plot_data, plot_colors = self._generate_plot_data(packets_mode)
        bin_edges = (self.velocity).to("km/s")

        if bins_range:
            bin_edges = bin_edges[bins_range[0] - 1 : bins_range[1] + 1]

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

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
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")
        if log_scale:
            self.ax.set_xscale("log")
            self.ax.set_yscale("log")
        plt.tight_layout()

        return self.ax

    def generate_plot_ply(
        self,
        packets_mode="virtual",
        fig=None,
        graph_height=500,
        cmapname="jet",
        species_list=None,
        log_scale=False,
        bins_range=None,
    ):
        """
        Generate the last interaction radius distribution plot using plotly.

        Parameters
        ----------
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        fig : plotly.graph_objects.Figure, optional
            Plotly figure object to add the plot to. If None, creates a new figure.
        graph_height : int, optional
            Height of the graph in pixels. Default is 500.
        cmapname : str, optional
            Colormap name. Default is 'jet'. A specific colormap can be chosen, such as "jet", "viridis", "plasma", etc.
        species_list : list of str
            List of species to plot.
        log_scale : bool, optional
            If True, both axes are scaled logarithmically. Default is False.
        bins_range : tuple of int, optional
            Range of bins to plot, e.g., (3, 9). Default is None, which plots all bins.

        Returns
        -------
        plotly.graph_objects.Figure
            Plotly figure object with the plot.
        """
        self._prepare_plot(packets_mode, species_list, cmapname)
        plot_data, plot_colors = self._generate_plot_data(packets_mode)
        bin_edges = (self.velocity).to("km/s")

        if bins_range:
            bin_edges = bin_edges[bins_range[0] - 1 : bins_range[1] + 1]

        fig = go.Figure()
        for data, color, name in zip(
            plot_data, plot_colors, self._species_name
        ):
            hist, _ = np.histogram(data, bins=bin_edges)
            step_x = np.repeat(bin_edges, 2)[1:-1]
            step_y = np.repeat(hist, 2)
            fig.add_trace(
                go.Scatter(
                    x=step_x,
                    y=step_y,
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
            height=graph_height,
            xaxis_title="Last Interaction Velocity (km/s)",
            yaxis_title="Packet Count",
            font=dict(size=14),
            yaxis=dict(exponentformat="power" if log_scale else "e"),
            xaxis=dict(exponentformat="power" if log_scale else "none"),
        )
        if log_scale:
            fig.update_xaxes(type="log")
            fig.update_yaxes(type="log", dtick=1)

        return fig
