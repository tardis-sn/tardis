import logging
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

logger = logging.getLogger(__name__)


class LIVPlotter:
    """
    Plotting interface for the last interaction velocity plot.
    """

    def __init__(self, data, time_explosion, velocity):
        """
        Initialize the plotter with required data from the simulation.

        Parameters
        ----------
        data : dict of SDECData
            Dictionary to store data required for last interaction velocity plot,
            for both packet modes (real, virtual).

        time_explosion : astropy.units.Quantity
            Time of the explosion.

        velocity : astropy.units.Quantity
            Velocity array from the simulation.
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
        LIVPlotter
        """

        return cls(
            dict(
                virtual=sdec.SDECData.from_simulation(sim, "virtual"),
                real=sdec.SDECData.from_simulation(sim, "real"),
            ),
            sim.plasma.time_explosion,
            sim.simulation_state.velocity,
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
        LIVPlotter
        """
        with pd.HDFStore(hdf_fpath, "r") as hdf:
            time_explosion = (
                hdf["/simulation/plasma/scalars"]["time_explosion"] * u.s
            )
            v_inner = hdf["/simulation/simulation_state/v_inner"] * (u.cm / u.s)
            v_outer = hdf["/simulation/simulation_state/v_outer"] * (u.cm / u.s)
            velocity = pd.concat(
                [v_inner, pd.Series([v_outer.iloc[-1]])], ignore_index=True
            ).tolist() * (u.cm / u.s)
            return cls(
                dict(
                    virtual=sdec.SDECData.from_hdf(hdf_fpath, "virtual"),
                    real=sdec.SDECData.from_hdf(hdf_fpath, "real"),
                ),
                time_explosion,
                velocity,
            )

    def _parse_species_list(self, species_list, packets_mode, nelements=None):
        """
        Parse user requested species list and create list of species ids to be used.

        Parameters
        ----------
        species_list : list of species to plot
            List of species (e.g. Si II, Ca II, etc.) that the user wants to show as unique colours.
            Species can be given as an ion (e.g. Si II), an element (e.g. Si), a range of ions
            (e.g. Si I - V), or any combination of these (e.g. species_list = [Si II, Fe I-V, Ca])
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        nelements : int, optional
            Number of elements to include in plot. The most interacting elements are included. If None, displays all elements.

        Raises
        ------
        ValueError
            If species list contains invalid entries.

        """
        self.sdec_plotter._parse_species_list(species_list)
        self._species_list = self.sdec_plotter._species_list
        self._species_mapped = self.sdec_plotter._species_mapped
        self._keep_colour = self.sdec_plotter._keep_colour

        if nelements:
            interaction_counts = (
                self.data[packets_mode]
                .packets_df_line_interaction["last_line_interaction_species"]
                .value_counts()
            )
            interaction_counts.index = interaction_counts.index // 100
            element_counts = interaction_counts.groupby(
                interaction_counts.index
            ).sum()
            top_elements = element_counts.nlargest(nelements).index
            top_species_list = [
                atomic_number2element_symbol(element)
                for element in top_elements
            ]
            self._parse_species_list(top_species_list, packets_mode)

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

        This method creates a list of colors corresponding to the species names.
        The colors are generated based on the species present in the model and
        the requested species list.
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
        """
        groups = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_species")
        )

        self.plot_colors = []
        self.plot_data = []
        species_not_wvl_range = []
        species_counter = 0

        for specie_list in self._species_mapped.values():
            full_v_last = []
            for specie in specie_list:
                if specie in self.species:
                    if specie not in groups.groups:
                        atomic_number = specie // 100
                        ion_number = specie % 100
                        ion_numeral = int_to_roman(ion_number + 1)
                        label = f"{atomic_number2element_symbol(atomic_number)} {ion_numeral}"
                        species_not_wvl_range.append(label)
                        continue
                    g_df = groups.get_group(specie)
                    r_last_interaction = (
                        g_df["last_interaction_in_r"].values * u.cm
                    )
                    v_last_interaction = (
                        r_last_interaction / self.time_explosion
                    ).to("km/s")
                    full_v_last.extend(v_last_interaction)
            if full_v_last:
                self.plot_data.append(full_v_last)
                self.plot_colors.append(self._color_list[species_counter])
                species_counter += 1

        if species_not_wvl_range:
            logger.info(
                "%s were not found in the provided wavelength range.",
                species_not_wvl_range,
            )

    def _prepare_plot_data(
        self,
        packets_mode,
        packet_wvl_range,
        species_list,
        cmapname,
        num_bins,
        nelements,
    ):
        """
        Prepare data and settings required for generating a plot.

        This method handles the common logic for preparing data and settings
        needed to generate both matplotlib and plotly plots. It parses the species
        list, generates color labels and colormap, and bins the velocity data.

        Parameters
        ----------
        packets_mode : str
            Packet mode, either 'virtual' or 'real'.
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA
        species_list : list of str
            List of species to plot. Species can be specified as an ion
            (e.g., Si II), an element (e.g., Si), a range of ions (e.g., Si I-V),
            or any combination of these.
        cmapname : str
            Name of the colormap to use. A specific colormap can be chosen, such
            as "jet", "viridis", "plasma", etc.
        num_bins : int, optional
            Number of bins for regrouping within the same range. If None,
            no regrouping is done.
        nelements : int, optional
            Number of elements to include in plot. The most interacting elements are included. If None, displays all elements.

        Raises
        ------
        ValueError
            If no species are provided for plotting, or if no valid species are
            found in the model.
        """
        if species_list is None:
            # Extract all unique elements from the packets data
            species_in_model = np.unique(
                self.data[packets_mode]
                .packets_df_line_interaction["last_line_interaction_species"]
                .values
            )
            species_list = [
                f"{atomic_number2element_symbol(specie // 100)}"
                for specie in species_in_model
            ]
        self._parse_species_list(species_list, packets_mode, nelements)
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

        if packet_wvl_range is None:
            self.packet_nu_line_range_mask = np.ones(
                self.data[packets_mode].packets_df_line_interaction.shape[0],
                dtype=bool,
            )
        else:
            packet_nu_range = [
                value.to("Hz", equivalencies=u.spectral())
                for value in packet_wvl_range
            ]
            self.packet_nu_line_range_mask = (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                >= packet_nu_range[1]
            ) & (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                <= packet_nu_range[0]
            )

        self._generate_plot_data(packets_mode)
        bin_edges = (self.velocity).to("km/s")

        if num_bins:
            if num_bins < 1:
                raise ValueError("Number of bins must be positive")
            elif num_bins > len(bin_edges) - 1:
                logger.warning(
                    "Number of bins must be less than or equal to number of shells. Plotting with number of bins equals to number of shells."
                )
                self.new_bin_edges = bin_edges
            else:
                self.new_bin_edges = np.linspace(
                    bin_edges[0], bin_edges[-1], num_bins + 1
                )
        else:
            self.new_bin_edges = bin_edges

    def _get_step_plot_data(self, data, bin_edges):
        """
        Generate step plot data from histogram data.

        Parameters
        ----------
        data : array-like
            Data to be binned into a histogram.
        bin_edges : array-like
            Edges of the bins for the histogram.
        """
        hist, _ = np.histogram(data, bins=bin_edges)
        self.step_x = np.repeat(bin_edges, 2)[1:-1]
        self.step_y = np.repeat(hist, 2)

    def generate_plot_mpl(
        self,
        species_list=None,
        nelements=None,
        packets_mode="virtual",
        packet_wvl_range=None,
        ax=None,
        figsize=(11, 5),
        cmapname="jet",
        xlog_scale=False,
        ylog_scale=False,
        num_bins=None,
        velocity_range=None,
    ):
        """
        Generate the last interaction velocity distribution plot using matplotlib.

        Parameters
        ----------
        species_list : list of str, optional
            List of species to plot. Default is None which plots all species in the model.
        nelements : int, optional
            Number of elements to include in plot. The most interacting elements are included. If None, displays all elements.
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, creates a new figure.
        figsize : tuple, optional
            Size of the figure. Default is (11, 5).
        cmapname : str, optional
            Colormap name. Default is 'jet'. A specific colormap can be chosen, such as "jet", "viridis", "plasma", etc.
        xlog_scale : bool, optional
            If True, x-axis is scaled logarithmically. Default is False.
        ylog_scale : bool, optional
            If True, y-axis is scaled logarithmically. Default is False.
        num_bins : int, optional
            Number of bins for regrouping within the same range. Default is None.
        velocity_range : tuple, optional
            Limits for the x-axis. If specified, overrides any automatically determined limits.

        Returns
        -------
        matplotlib.axes.Axes
            Axes object with the plot.
        """
        # If species_list and nelements requested, tell user that nelements is ignored
        if species_list is not None and nelements is not None:
            logger.info(
                "Both nelements and species_list were requested. Species_list takes priority; nelements is ignored"
            )
            nelements = None

        self._prepare_plot_data(
            packets_mode,
            packet_wvl_range,
            species_list,
            cmapname,
            num_bins,
            nelements,
        )

        bin_edges = self.new_bin_edges

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        for data, color, name in zip(
            self.plot_data, self.plot_colors, self._species_name
        ):
            self._get_step_plot_data(data, bin_edges)
            self.ax.plot(
                self.step_x,
                self.step_y,
                label=name,
                color=color,
                linewidth=2.5,
                drawstyle="steps-post",
                alpha=0.75,
            )

        self.ax.ticklabel_format(axis="y", scilimits=(0, 0))
        self.ax.tick_params("both", labelsize=15)
        self.ax.set_xlabel("Last Interaction Velocity (km/s)", fontsize=14)
        self.ax.set_ylabel("Packet Count", fontsize=15)
        self.ax.legend(fontsize=15, bbox_to_anchor=(1.0, 1.0), loc="upper left")
        self.ax.figure.tight_layout()
        if xlog_scale:
            self.ax.set_xscale("log")
        if ylog_scale:
            self.ax.set_yscale("log")
        if velocity_range:
            self.ax.set_xlim(velocity_range[0], velocity_range[1])

        return self.ax

    def generate_plot_ply(
        self,
        species_list=None,
        nelements=None,
        packets_mode="virtual",
        packet_wvl_range=None,
        fig=None,
        graph_height=600,
        cmapname="jet",
        xlog_scale=False,
        ylog_scale=False,
        num_bins=None,
        velocity_range=None,
    ):
        """
        Generate the last interaction velocity distribution plot using plotly.

        Parameters
        ----------
        species_list : list of str, optional
            List of species to plot. Default is None which plots all species in the model.
        nelements : int, optional
            Number of elements to include in plot. The most interacting elements are included. If None, displays all elements.
        packets_mode : str, optional
            Packet mode, either 'virtual' or 'real'. Default is 'virtual'.
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA
        fig : plotly.graph_objects.Figure, optional
            Plotly figure object to add the plot to. If None, creates a new figure.
        graph_height : int, optional
            Height (in px) of the plotly graph to display. Default value is 600.
        cmapname : str, optional
            Colormap name. Default is 'jet'. A specific colormap can be chosen, such as "jet", "viridis", "plasma", etc.
        xlog_scale : bool, optional
            If True, x-axis is scaled logarithmically. Default is False.
        ylog_scale : bool, optional
            If True, y-axis is scaled logarithmically. Default is False.
        num_bins : int, optional
            Number of bins for regrouping within the same range. Default is None.
        velocity_range : tuple, optional
            Limits for the x-axis. If specified, overrides any automatically determined limits.

        Returns
        -------
        plotly.graph_objects.Figure
            Plotly figure object with the plot.
        """
        # If species_list and nelements requested, tell user that nelements is ignored
        if species_list is not None and nelements is not None:
            logger.info(
                "Both nelements and species_list were requested. Species_list takes priority; nelements is ignored"
            )
            nelements = None

        self._prepare_plot_data(
            packets_mode,
            packet_wvl_range,
            species_list,
            cmapname,
            num_bins,
            nelements,
        )

        bin_edges = self.new_bin_edges

        if fig is None:
            self.fig = go.Figure()
        else:
            self.fig = fig

        for data, color, name in zip(
            self.plot_data, self.plot_colors, self._species_name
        ):
            self._get_step_plot_data(data, bin_edges)
            self.fig.add_trace(
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
                )
            )
        self.fig.update_layout(
            height=graph_height,
            xaxis_title="Last Interaction Velocity (km/s)",
            yaxis_title="Packet Count",
            font=dict(size=15),
            yaxis=dict(exponentformat="power" if ylog_scale else "e"),
            xaxis=dict(exponentformat="power" if xlog_scale else "none"),
        )
        if xlog_scale:
            self.fig.update_xaxes(type="log")
        if ylog_scale:
            self.fig.update_yaxes(type="log", dtick=1)

        if velocity_range:
            self.fig.update_xaxes(range=velocity_range)

        return self.fig
