import tardis.visualization.tools.sdec_plot as sdec

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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go


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
        self.no_of_shells = no_of_shells
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
        Plotter
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
        Plotter
        """

        return cls(
            dict(
                virtual=sdec.SDECData.from_hdf(hdf_fpath, "virtual"),
                real=sdec.SDECData.from_hdf(hdf_fpath, "real"),
            ),
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

        """
        if species_list is not None:
            # check if there are any digits in the species list. If there are, then exit.
            # species_list should only contain species in the Roman numeral
            # format, e.g. Si II, and each ion must contain a space
            if any(char.isdigit() for char in " ".join(species_list)) == True:
                raise ValueError(
                    "All species must be in Roman numeral form, e.g. Si II"
                )
            else:
                full_species_list = []
                for species in species_list:
                    # check if a hyphen is present. If it is, then it indicates a
                    # range of ions. Add each ion in that range to the list as a new entry
                    if "-" in species:
                        # split the string on spaces. First thing in the list is then the element
                        element = species.split(" ")[0]
                        # Next thing is the ion range
                        # convert the requested ions into numerals
                        first_ion_numeral = roman_to_int(
                            species.split(" ")[-1].split("-")[0]
                        )
                        second_ion_numeral = roman_to_int(
                            species.split(" ")[-1].split("-")[-1]
                        )
                        # add each ion between the two requested into the species list
                        for ion_number in np.arange(
                            first_ion_numeral, second_ion_numeral + 1
                        ):
                            full_species_list.append(
                                f"{element} {int_to_roman(ion_number)}"
                            )
                    else:
                        # Otherwise it's either an element or ion so just add to the list
                        full_species_list.append(species)

                # full_species_list is now a list containing each individual species requested
                # e.g. it parses species_list = [Si I - V] into species_list = [Si I, Si II, Si III, Si IV, Si V]
                self._full_species_list = full_species_list
                requested_species_ids = []
                keep_colour = []

                # go through each of the requested species. Check whether it is
                # an element or ion (ions have spaces). If it is an element,
                # add all possible ions to the ions list. Otherwise just add
                # the requested ion
                for species in full_species_list:
                    if " " in species:
                        requested_species_ids.append(
                            [
                                species_string_to_tuple(species)[0] * 100
                                + species_string_to_tuple(species)[1]
                            ]
                        )
                    else:
                        atomic_number = element_symbol2atomic_number(species)
                        requested_species_ids.append(
                            [
                                atomic_number * 100 + ion_number
                                for ion_number in np.arange(atomic_number)
                            ]
                        )
                        # add the atomic number to a list so you know that this element should
                        # have all species in the same colour, i.e. it was requested like
                        # species_list = [Si]
                        keep_colour.append(atomic_number)
                requested_species_ids = [
                    species_id
                    for temp_list in requested_species_ids
                    for species_id in temp_list
                ]

                self._species_list = requested_species_ids
                self._keep_colour = keep_colour
        else:
            self._species_list = None
        return

    def _make_colorbar_labels(self):
        """Get the labels for the species in the colorbar."""
        if self._species_list is None:
            # If species_list is none then the labels are just elements
            species_name = [
                atomic_number2element_symbol(atomic_num)
                for atomic_num in self.species
            ]
        else:
            species_name = []
            for species in self.species:
                # Go through each species requested
                ion_number = species % 100
                atomic_number = (species - ion_number) / 100

                ion_numeral = int_to_roman(ion_number + 1)
                atomic_symbol = atomic_number2element_symbol(atomic_number)

                # if the element was requested, and not a specific ion, then
                # add the element symbol to the label list
                if (atomic_number in self._keep_colour) & (
                    atomic_symbol not in species_name
                ):
                    # compiling the label, and adding it to the list
                    label = f"{atomic_symbol}"
                    species_name.append(label)
                elif atomic_number not in self._keep_colour:
                    # otherwise add the ion to the label list
                    label = f"{atomic_symbol} {ion_numeral}"
                    species_name.append(label)

        self._species_name = species_name
        return

    def _make_colorbar_colors(self):
        """Get the colours for the species to be plotted."""
        # the colours depends on the species present in the model and what's requested
        # some species need to be shown in the same colour, so the exact colours have to be
        # worked out

        color_list = []

        # Colors for each element
        # Create new variables to keep track of the last atomic number that was plotted
        # This is used when plotting species in case an element was given in the list
        # This is to ensure that all ions of that element are grouped together
        # ii is to track the colour index
        # e.g. if Si is given in species_list, this is to ensure Si I, Si II, etc. all have the same colour
        color_counter = 0
        previous_atomic_number = 0
        for species_counter, identifier in enumerate(self.species):
            if self._species_list is not None:
                # Get the ion number and atomic number for each species
                ion_number = identifier % 100
                atomic_number = (identifier - ion_number) / 100
                if previous_atomic_number == 0:
                    # If this is the first species being plotted, then take note of the atomic number
                    # don't update the colour index
                    color_counter = color_counter
                    previous_atomic_number = atomic_number
                elif previous_atomic_number in self._keep_colour:
                    # If the atomic number is in the list of elements that should all be plotted in the same colour
                    # then don't update the colour index if this element has been plotted already
                    if previous_atomic_number == atomic_number:
                        color_counter = color_counter
                        previous_atomic_number = atomic_number
                    else:
                        # Otherwise, increase the colour counter by one, because this is a new element
                        color_counter = color_counter + 1
                        previous_atomic_number = atomic_number
                else:
                    # If this is just a normal species that was requested then increment the colour index
                    color_counter = color_counter + 1
                    previous_atomic_number = atomic_number
                # Calculate the colour of this species
                color = self.cmap(color_counter / len(self._species_name))

            else:
                # If you're not using species list then this is just a fraction based on the total
                # number of columns in the dataframe
                color = self.cmap(species_counter / len(self.species))

            color_list.append(color)

        self._color_list = color_list

        return

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
        figsize=(12, 7),
        cmapname="jet",
        species_list=None,
    ):
        """
        Generate the last interaction radius distribution plot
        using matplotlib.
        """

        # Parse the requested species list
        self._parse_species_list(species_list=species_list)
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
        self._show_colorbar_mpl()

        groups = self.data[packets_mode].packets_df_line_interaction.groupby(
            by="last_line_interaction_species"
        )

        plot_colors = []
        plot_data = []

        for species_counter, identifier in enumerate(self.species):
            g_df = groups.get_group(identifier)
            r_last_interaction = g_df["last_interaction_in_r"].values * u.cm
            v_last_interaction = (r_last_interaction / self.time_explosion).to(
                "km/s"
            )
            plot_data.append(v_last_interaction)
            plot_colors.append(self._color_list[species_counter])

        self.ax.hist(plot_data, bins=self.no_of_shells, color=plot_colors)
        self.ax.ticklabel_format(axis="y", scilimits=(0, 0))
        self.ax.tick_params("both", labelsize=20)
        self.ax.set_xlabel("Last Interaction Velocity (km/s)", fontsize=25)
        self.ax.set_ylabel("Packet Count", fontsize=25)

        return plt.gca()

    def generate_plot_plotly(
        self,
        packets_mode="virtual",
        species_list=None,
    ):
        """
        Generate the last interaction radius distribution plot
        using Plotly.
        """
        # Parse the requested species list
        self._parse_species_list(species_list=species_list)
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
        # Get the number of unique colors
        self._make_colorbar_colors()
        groups = self.data[packets_mode].packets_df_line_interaction.groupby(
            by="last_line_interaction_species"
        )

        plot_colors = []
        plot_data = []

        for species_counter, identifier in enumerate(self.species):
            g_df = groups.get_group(identifier)
            r_last_interaction = g_df["last_interaction_in_r"].values * u.cm
            v_last_interaction = (r_last_interaction / self.time_explosion).to(
                "km/s"
            )
            plot_data.append(v_last_interaction)
            color = f"rgba({int(255*self._color_list[species_counter][0])}, {int(255*self._color_list[species_counter][1])}, {int(255*self._color_list[species_counter][2])}, 1)"
            plot_colors.append(color)
        fig = go.Figure()
        for data, color, name in zip(
            plot_data, plot_colors, self._species_name
        ):
            fig.add_trace(
                go.Histogram(
                    x=data,
                    marker_color=color,
                    name=name,
                    opacity=0.75,
                    nbinsx=self.no_of_shells,
                )
            )
        fig.update_layout(
            barmode="overlay",
            xaxis_title="Last Interaction Velocity (km/s)",
            yaxis_title="Packet Count",
            template="plotly_white",
            font=dict(size=18),
        )
        return fig
