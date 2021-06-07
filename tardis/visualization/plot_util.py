"""Utility functions to be used in plotting."""

import re
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go

from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    species_string_to_tuple,
    roman_to_int,
    int_to_roman,
)


def axis_label_in_latex(label_text, unit, only_text=True):
    """
    Get axis label for plotly plots that can show units in latex.

    Parameters
    ----------
    label_text : str
        Text to show on label, may be expressed in latex
    unit : astropy.units
        Unit of the label which needs to be expressed in latex
    only_text : bool
        If label_text is expressed purely in text (i.e. without
        using latex) or not. Default value is True

    Returns
    -------
    str
        Latex string for label renderable by plotly
    """
    unit_in_latex = unit.to_string("latex_inline").strip("$")

    # If present, place s^{-1} just after erg
    if "erg" in unit_in_latex and "s^{-1}" in unit_in_latex:
        constituent_units = (
            re.compile("\\\mathrm\{(.*)\}")
            .findall(unit_in_latex)[0]
            .split("\\,")
        )
        constituent_units.remove("s^{-1}")
        constituent_units.insert(constituent_units.index("erg") + 1, "s^{-1}")
        constituent_units_string = "\\,".join(constituent_units)
        unit_in_latex = f"\\mathrm{{{constituent_units_string}}}"

    if only_text:
        return f"$\\text{{{label_text}}}\\,[{unit_in_latex}]$"
    else:
        return f"${label_text}\\,[{unit_in_latex}]$"


def get_mid_point_idx(arr):
    """
    Get index of the middle point of a sorted array (ascending or descending).

    The values in array may not be evenly distributed so it picks the middle
    point not by index but by their values.

    Parameters
    ----------
    arr : np.array

    Returns
    -------
    int
    """
    mid_value = (arr[0] + arr[-1]) / 2
    return np.abs(arr - mid_value).argmin()


class Plotter:
    """
    Class to encapsulate common functions used by other plotter tools.
    """

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

    def _show_colorbar_ply(self):
        """Show plotly colorbar with labels of elements mapped to colors."""
        # Interpolate [0, 1] range to create bins equal to number of elements
        colorscale_bins = np.linspace(0, 1, num=len(self._species_name) + 1)

        # Create a categorical colorscale [a list of (reference point, color)]
        # by mapping same reference points (excluding 1st and last bin edge)
        # twice in a row (https://plotly.com/python/colorscales/#constructing-a-discrete-or-discontinuous-color-scale)
        categorical_colorscale = []
        for species_counter in range(len(self._species_name)):
            color = self.to_rgb255_string(
                self.cmap(colorscale_bins[species_counter])
            )
            categorical_colorscale.append(
                (colorscale_bins[species_counter], color)
            )
            categorical_colorscale.append(
                (colorscale_bins[species_counter + 1], color)
            )

        coloraxis_options = dict(
            colorscale=categorical_colorscale,
            showscale=True,
            cmin=0,
            cmax=len(self._species_name),
            colorbar=dict(
                title="Elements",
                tickvals=np.arange(0, len(self._species_name)) + 0.5,
                ticktext=self._species_name,
                # to change length and position of colorbar
                len=0.75,
                yanchor="top",
                y=0.75,
            ),
        )

        # Plot an invisible one point scatter trace, to make colorbar show up
        scatter_point_idx = get_mid_point_idx(self.plot_wavelength)
        self.fig.add_trace(
            go.Scatter(
                x=self.plot_wavelength[scatter_point_idx],
                y=[0],
                mode="markers",
                showlegend=False,
                hoverinfo="skip",
                marker=dict(color=[0], opacity=0, **coloraxis_options),
            )
        )
