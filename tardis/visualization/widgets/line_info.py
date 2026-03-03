"""Class to create and display Line Info Widget."""

import logging

import numpy as np
import pandas as pd
import panel as pn
import matplotlib.pyplot as plt
from astropy import units as u
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure

from tardis.analysis import LastLineInteraction
from tardis.configuration.sorting_globals import SORTING_ALGORITHM
from tardis.util.base import (
    atomic_number2element_symbol,
    species_string_to_tuple,
    species_tuple_to_string,
)
from tardis.util.environment import Environment
from tardis.visualization import plot_util as pu
from tardis.visualization.tools.sdec_plot import SDECPlotter
from tardis.visualization.widgets.util import (
    TableSummaryLabel,
    create_table_widget,
)

logger = logging.getLogger(__name__)

class LineInfoWidget:
    """
    Widget to explore atomic lines that produced features in the simulated spectrum.

    It allows selection of a wavelength range in the spectrum to display a
    table giving the fraction of packets that experienced their last
    interaction with each species. Using toggle buttons, users can specify
    whether to filter the selected range by emitted or absorbed wavelengths
    of packets. Clicking on a row in the fractional species interactions table
    shows packet counts for each last line interaction of the selected species,
    which can be grouped in several ways using the dropdown menu.
    """

    FILTER_MODES = ("packet_out_nu", "packet_in_nu")
    FILTER_MODES_DESC = ("Emitted Wavelength", "Absorbed Wavelength")
    GROUP_MODES = ("both", "exc", "de-exc")
    GROUP_MODES_DESC = (
        "Both excitation line (absorption) and de-excitation line (emission)",
        "Only excitation line (absorption)",
        "Only de-excitation line (emission)",
    )
    COLORS = {"selection_area": "lightpink", "selection_border": "salmon"}

    def __init__(
            self,
            lines_data,
            line_interaction_analysis,
            spectrum_wavelength,
            spectrum_luminosity_density_lambda,
            virt_spectrum_wavelength,
            virt_spectrum_luminosity_density_lambda,
            sdec_plotter=None,
            show_sdec=False,
            sdec_packets_mode="virtual",
            observed_spectrum=None,
    ):
        """
        Initialize the LineInfoWidget with line interaction and spectrum data.

        Parameters
        ----------
        lines_data : pd.DataFrame
            Data about the atomic lines present in simulation model's plasma
        line_interaction_analysis : dict of tardis.analysis.LastLineInteraction
            Dictionary in which keys are the FILTER_MODES and values are the
            LastLineInteraction objects initialized with corresponding modes
        spectrum_wavelength : astropy.Quantity
            Wavelength values of a real spectrum, having unit of Angstrom
        spectrum_luminosity_density_lambda : astropy.Quantity
            Luminosity density lambda values of a real spectrum, having unit
            of (erg/s)/Angstrom
        virt_spectrum_wavelength : astropy.Quantity
            Wavelength values of a virtual spectrum, having unit of Angstrom
        virt_spectrum_luminosity_density_lambda : astropy.Quantity
            Luminosity density lambda values of a virtual spectrum, having unit
            of (erg/s)/Angstrom
        sdec_plotter : SDECPlotter, optional
            SDECPlotter instance for displaying SDEC plot data
        show_sdec : bool, optional
            Whether to plot SDEC data (default: False)
        sdec_packets_mode : str, optional
            Packets mode to use for SDEC plotter (default: "virtual")
        observed_spectrum : tuple or list of astropy.Quantity, optional
            Option to plot an observed spectrum in the widget. If given, the first element
            should be the wavelength and the second element should be flux,
            i.e. (wavelength, flux). The assumed units for wavelength and flux are
            angstroms and erg/(angstroms * s * cm^2), respectively. Default value is None.
        """
        self.lines_data = lines_data
        self.line_interaction_analysis = line_interaction_analysis
        self.sdec_plotter = sdec_plotter
        self.show_sdec = show_sdec
        self.sdec_packets_mode = sdec_packets_mode
        self.observed_spectrum = observed_spectrum

        # Store renderers for toggle functionality
        self.line_renderers = {}

        self.element_checkboxes = {}
        self.checkbox_real_packets = pn.widgets.Checkbox(name="Real Packets", value=True)
        self.checkbox_virtual_packets = pn.widgets.Checkbox(name="Virtual Packets", value=True)
        if self.show_sdec:
            self.checkbox_sdec_spectrum = pn.widgets.Checkbox(name="SDEC Spectrum", value=True)
            self.checkbox_photosphere = pn.widgets.Checkbox(name="Blackbody Photosphere", value=True)
            self.checkbox_no_interaction = pn.widgets.Checkbox(name="No Interaction", value=True)
            self.checkbox_electron_scatter = pn.widgets.Checkbox(name="Electron Scatter Only", value=True)
        if self.observed_spectrum is not None:
            self.checkbox_observed_spectrum = pn.widgets.Checkbox(name="Observed Spectrum", value=True)

        # Link to callbacks
        checkboxes = [self.checkbox_real_packets, self.checkbox_virtual_packets]
        if self.show_sdec:
            checkboxes.extend([
                self.checkbox_sdec_spectrum, self.checkbox_photosphere,
                self.checkbox_no_interaction, self.checkbox_electron_scatter
            ])
        if self.observed_spectrum is not None:
            checkboxes.append(self.checkbox_observed_spectrum)
        for cb in checkboxes:
            cb.param.watch(self._on_toggle_line, "value")

        # Widgets ------------------------------------------------
        max_rows_option = {"maxVisibleRows": 9}
        self.species_interactions_table = create_table_widget(
            data=self.get_species_interactions(None),
            table_options=max_rows_option,
        )

        self.last_line_counts_table = create_table_widget(
            data=self.get_last_line_counts(None),
            table_options=max_rows_option,
        )
        self.total_packets_label = TableSummaryLabel(
            target_table=self.last_line_counts_table,
            table_col_widths=[75, 25],
            label_key="Total Packets",
            label_value=0,
        )

        self.figure_widget = self.plot_spectrum(
            spectrum_wavelength,
            spectrum_luminosity_density_lambda,
            virt_spectrum_wavelength,
            virt_spectrum_luminosity_density_lambda,
        )

        self.filter_mode_buttons = pn.widgets.RadioButtonGroup(
            options=list(self.FILTER_MODES_DESC),
            value=self.FILTER_MODES_DESC[0],
        )

        self.group_mode_dropdown = pn.widgets.Select(
            options=list(self.GROUP_MODES_DESC), value=self.GROUP_MODES_DESC[0]
        )

        self._current_wavelength_range = None  # Track current selection

    @classmethod
    def from_simulation(cls, sim, show_sdec=False, sdec_packets_mode="virtual",  observed_spectrum=None):
        """
        Create an instance of LineInfoWidget from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
        show_sdec : bool, optional
            Whether to plot SDEC data (default: False)
        observed_spectrum : tuple or list of astropy.Quantity, optional
            Option to plot an observed spectrum in the widget. If given, the first element
            should be the wavelength and the second element should be flux,
            i.e. (wavelength, flux). The assumed units for wavelength and flux are
            angstroms and erg/(angstroms * s * cm^2), respectively. Default value is None.

        Returns
        -------
        LineInfoWidget object
        """
        spectrum_solver = sim.spectrum_solver
        sdec_plotter = SDECPlotter.from_simulation(sim) if show_sdec else None

        return cls(
            lines_data=sim.plasma.lines.reset_index().set_index("line_id"),
            line_interaction_analysis={
                filter_mode: LastLineInteraction.from_simulation(
                    sim, filter_mode
                )
                for filter_mode in cls.FILTER_MODES
            },
            spectrum_wavelength=spectrum_solver.spectrum_real_packets.wavelength,
            spectrum_luminosity_density_lambda=spectrum_solver.spectrum_real_packets.luminosity_density_lambda.to(
                "erg/(s AA)"
            ),
            virt_spectrum_wavelength=spectrum_solver.spectrum_virtual_packets.wavelength,
            virt_spectrum_luminosity_density_lambda=spectrum_solver.spectrum_virtual_packets.luminosity_density_lambda.to(
                "erg/(s AA)"
            ),
            sdec_plotter=sdec_plotter,
            show_sdec=show_sdec,
            sdec_packets_mode=sdec_packets_mode,
            observed_spectrum=observed_spectrum,
        )

    def get_species_interactions(
        self, wavelength_range, filter_mode=FILTER_MODES[0]
    ):
        """
        Get fractional species interactions in specified wavelength range.

        Fractional species interactions means fraction of packets present in
        the specified wavelength range which experienced their last interaction
        with a species. The packets to consider are filtered by the specified
        filter mode.

        Parameters
        ----------
        wavelength_range : list-like or None
            A list of two float values to specify the wavelength range - first
            for the range start and second for the range end. None specifies
            that no wavelength range is selected and will return empty dataframe
        filter_mode : str, optional
            Filter mode of the LastLineInteraction object to use for filtering
            the selected wavelength range (more details in Notes section).
            Allowed values are given by the class variable :code:`FILTER_MODES`
            (default value is :code:`FILTER_MODES[0]`)

        Returns
        -------
        pandas.DataFrame
            Dataframe containing species symbols and corresponding fractions
            of packets interacting with them

        Notes
        -----
        This method depends on tardis.analysis.LastLineInteraction object for
        doing computations. So there is a member variable in this class -
        :code:`line_interaction_analysis` which is a dictionary of such objects
        (each of them differ in how they filter the selected wavelength range).
        Thus we have to specify which object to use by specifying the
        filter_mode parameter.
        """
        if wavelength_range:
            self.line_interaction_analysis[filter_mode].wavelength_start = (
                wavelength_range[0] * u.AA
            )
            self.line_interaction_analysis[filter_mode].wavelength_end = (
                wavelength_range[1] * u.AA
            )

            # Obtain species group from last_line_in dataframe
            selected_species_group = self.line_interaction_analysis[
                filter_mode
            ].last_line_in.groupby(["atomic_number", "ion_number"])

            if selected_species_group.groups:
                selected_species_symbols = [
                    species_tuple_to_string(item)
                    for item in selected_species_group.groups.keys()
                ]

                # Normalize each group's count to find fractions of interactions
                fractional_species_interactions = (
                    selected_species_group.size()
                    / self.line_interaction_analysis[
                        filter_mode
                    ].last_line_in.shape[0]
                )

            else:  # No species could be selected in specified wavelength_range
                # Create one row with empty strings for empty dataframe
                selected_species_symbols = [""]
                fractional_species_interactions = pd.Series([""])

        else:  # wavelength_range is None
            selected_species_symbols = [""]
            fractional_species_interactions = pd.Series([""])

        fractional_species_interactions.index = pd.Index(
            selected_species_symbols, name="Species"
        )
        fractional_species_interactions.name = "Fraction of packets interacting"
        return fractional_species_interactions.sort_values(
            ascending=False, kind=SORTING_ALGORITHM
        ).to_frame()

    def get_last_line_counts(
        self,
        selected_species,
        filter_mode=FILTER_MODES[0],
        group_mode=GROUP_MODES[0],
    ):
        """
        Get packet counts of each last line interaction of a species.

        Parameters
        ----------
        selected_species : str
            Valid symbol of a species (e.g Si II) selected from the species
            data returned by :code:`get_species_interactions` (see Notes section)
        filter_mode : str, optional
            Filter mode of the LastLineInteraction object to use for fetching
            the data of last lines interacted (more details in Notes section).
            Allowed values are given by the class variable :code:`FILTER_MODES`
            (default value is :code:`FILTER_MODES[0]`)
        group_mode : str, optional
            Group mode to use for grouping last line interactions by excitation
            lines, de-excitation lines or both. Allowed values are given by the
            class variable :code:`GROUP_MODES` (default value is
            :code:`GROUP_MODES[0]` i.e. both)

        Returns
        -------
        pd.DataFrame
            DataFrame containing last line interactions and corresponding
            packet counts.

        Notes
        -----
        This method depends on tardis.analysis.LastLineInteraction object for
        doing computations. So there is a member variable in this class -
        :code:`line_interaction_analysis` which is a dictionary of such objects
        (each of them differ in how they filter the selected wavelength range).
        Thus we have to specify which object to use by specifying the
        filter_mode parameter.

        This method should always be called after calling
        :code:`get_species_interactions` method which sets a wavelength
        range on LastLineInteraction object. So selected_species should
        be one present within that range, otherwise it will result an error.
        """
        if selected_species:
            selected_species_tuple = species_string_to_tuple(selected_species)

            try:
                # Get selected species' rows from last_line_in dataframe
                current_last_lines_in = (
                    self.line_interaction_analysis[filter_mode]
                    .last_line_in.xs(
                        key=selected_species_tuple,
                        level=["atomic_number", "ion_number"],
                        drop_level=False,
                    )
                    .reset_index()
                )

                # Get selected species' rows from last_line_out dataframe
                current_last_lines_out = (
                    self.line_interaction_analysis[filter_mode]
                    .last_line_out.xs(
                        key=selected_species_tuple,
                        level=["atomic_number", "ion_number"],
                        drop_level=False,
                    )
                    .reset_index()
                )

                assert (
                    current_last_lines_in.empty & current_last_lines_out.empty
                    is False
                )

            except (KeyError, AssertionError):  # selected_species is invalid
                allowed_species = [
                    species_tuple_to_string(species)
                    for species in self.line_interaction_analysis[filter_mode]
                    .last_line_in.groupby(["atomic_number", "ion_number"])
                    .groups.keys()
                ]
                raise ValueError(
                    "Invalid value of selected_species, it must be one present "
                    "within the currently selected wavelength range in your "
                    f"LineInfoWidget instance, which are {allowed_species}"
                )

            last_line_interaction_string = []
            interacting_packets_count = []

            if group_mode == "both":
                # Group by both exc. line ids and de-exc. line ids
                current_last_lines_in["line_id_out"] = (
                    current_last_lines_out.line_id
                )
                grouped_line_interactions = current_last_lines_in.groupby(
                    ["line_id", "line_id_out"]
                )

                # Iterate over each group's key and size and append them to list
                for (
                    line_id,
                    count,
                ) in grouped_line_interactions.size().items():
                    current_line_in = self.lines_data.loc[line_id[0]]
                    current_line_out = self.lines_data.loc[line_id[1]]
                    last_line_interaction_string.append(
                        f"exc. {int(current_line_in.level_number_lower):02d}-"
                        f"{int(current_line_in.level_number_upper):02d} "
                        f"({current_line_in.wavelength:.2f} A) "
                        f"de-exc. {int(current_line_out.level_number_upper):02d}-"
                        f"{int(current_line_out.level_number_lower):02d} "
                        f"({current_line_out.wavelength:.2f} A)"
                    )
                    interacting_packets_count.append(count)

            elif group_mode == "exc":
                grouped_line_interactions = current_last_lines_in.groupby(
                    "line_id"
                )

                # Iterate over each group's key and size and append them to list
                for (
                    line_id,
                    count,
                ) in grouped_line_interactions.size().items():
                    current_line_in = self.lines_data.loc[line_id]
                    last_line_interaction_string.append(
                        f"exc. {int(current_line_in.level_number_lower):02d}-"
                        f"{int(current_line_in.level_number_upper):02d} "
                        f"({current_line_in.wavelength:.2f} A)"
                    )
                    interacting_packets_count.append(count)

            elif group_mode == "de-exc":
                grouped_line_interactions = current_last_lines_out.groupby(
                    "line_id"
                )

                # Iterate over each group's key and size and append them to list
                for (
                    line_id,
                    count,
                ) in grouped_line_interactions.size().items():
                    current_line_out = self.lines_data.loc[line_id]
                    last_line_interaction_string.append(
                        f"de-exc. {int(current_line_out.level_number_upper):02d}-"
                        f"{int(current_line_out.level_number_lower):02d} "
                        f"({current_line_out.wavelength:.2f} A)"
                    )
                    interacting_packets_count.append(count)

            else:
                raise ValueError(
                    "Invalid value passed to group_mode argument. "
                    f"Allowed values are {self.GROUP_MODES}"
                )

        else:  # species_selected is None
            # Create one row with empty strings for empty dataframe
            interacting_packets_count = [""]
            last_line_interaction_string = [""]

        last_line_counts = pd.Series(interacting_packets_count)
        last_line_counts.name = "No. of packets"
        last_line_counts.index = pd.Index(
            last_line_interaction_string, name="Last Line Interaction"
        )
        return last_line_counts.sort_values(
            ascending=False, kind=SORTING_ALGORITHM
        ).to_frame()

    @staticmethod
    def get_middle_half_edges(arr):
        """
        Get edges of the middle half range of an array.

        Parameters
        ----------
        arr : np.array

        Returns
        -------
        list
        """
        arr = np.sort(arr, kind=SORTING_ALGORITHM)
        return [
            (arr[-1] - arr[0]) / 4 + arr[1],
            (arr[-1] - arr[0]) * 3 / 4 + arr[1],
        ]

    def plot_spectrum(
            self,
            wavelength,
            luminosity_density_lambda,
            virt_wavelength,
            virt_luminosity_density_lambda,
    ):
        """
        Produce a Bokeh figure by plotting the spectrum of model with SDEC data.
        """
        p = figure(
            width=800, height=400, title="Spectrum", tools="box_select,reset"
        )

        p.xaxis.axis_label = f"Wavelength [{wavelength.unit}]"
        p.yaxis.axis_label = f"Luminosity [{luminosity_density_lambda.unit}]"

        self.line_renderers["real_packets"] = p.line(
            wavelength.value,
            luminosity_density_lambda.value,
            legend_label="Real packets",
            color="blue",
        )

        self.line_renderers["virtual_packets"] = p.line(
            virt_wavelength.value,
            virt_luminosity_density_lambda.value,
            legend_label="Virtual packets",
            color="red",
        )

        if self.observed_spectrum is not None:
            observed_wavelength = self.observed_spectrum[0].to(u.AA)
            observed_flux = self.observed_spectrum[1]
            self.line_renderers["observed"] = p.line(
                observed_wavelength.value,
                observed_flux.value,
                legend_label="Observed Spectrum",
                color="black",
            )

        if self.show_sdec and self.sdec_plotter is not None:
            self.sdec_plotter.prepare_plot_data(
                packets_mode=self.sdec_packets_mode,
                packet_wvl_range=None,
                distance=None,
                species_list=None,
                nelements=None,
            )

            self.line_renderers["sdec_spectrum"] = p.line(
                self.sdec_plotter.plot_wavelength.value,
                self.sdec_plotter.modeled_spectrum_luminosity.value,
                legend_label=f"{self.sdec_packets_mode.capitalize()} Spectrum",
                color="green",
                line_dash="dashed",
            )

            self.line_renderers["photosphere"] = p.line(
                self.sdec_plotter.plot_wavelength.value,
                self.sdec_plotter.photosphere_luminosity.value,
                legend_label="Blackbody Photosphere",
                color="orange",
                line_dash="dotted",
            )

            self.line_renderers["no_interaction"] = p.line(
                self.sdec_plotter.plot_wavelength.value,
                self.sdec_plotter.emission_luminosities_df[("noint", "")].values,
                legend_label="No Interaction",
                color="#4C4C4C",
                line_width=1.5,
            )

            self.line_renderers["electron_scatter"] = p.line(
                self.sdec_plotter.plot_wavelength.value,
                self.sdec_plotter.emission_luminosities_df[("escatter", "")].values,
                legend_label="Electron Scatter Only",
                color="#8F8F8F",
                line_width=1.5,
            )

            self._plot_species_contributions(p)

        # Create invisible scatter for selection
        source = ColumnDataSource(
            dict(x=wavelength.value, y=luminosity_density_lambda.value)
        )
        p.scatter("x", "y", source=source, alpha=0, size=1)

        # Create selection overlay source (initially empty)
        selection_source = ColumnDataSource(
            dict(left=[], right=[], top=[], bottom=[])
        )
        p.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=selection_source,
            alpha=0.3,
            color="lightblue",
        )

        # Store references for callback
        self._bokeh_plot = p
        self._selection_source = selection_source
        self._y_range = [
            luminosity_density_lambda.value.min(),
            luminosity_density_lambda.value.max(),
        ]
        self._wavelength_data = (
            wavelength  # Store wavelength data for callback access
        )

        # Connect selection callback
        source.selected.on_change("indices", self._selection_callback)

        return pn.pane.Bokeh(p)

    def _plot_species_contributions(self, p):
        """Plot emission and absorption contributions for each element as stacked areas."""
        if not self.show_sdec or self.sdec_plotter is None:
            return

        emission_df = self.sdec_plotter.emission_luminosities_df
        absorption_df = self.sdec_plotter.absorption_luminosities_df
        wavelength = self.sdec_plotter.plot_wavelength.value
        species = self.sdec_plotter.species

        if len(species) == 0:
            return

        cmap = plt.get_cmap("jet", len(species))
        color_list = [
            pu.to_rgb255_string(cmap(i / len(species)))
            for i in range(len(species))
        ]
        self._element_info_for_legend = []

        emission_lower = np.zeros(len(wavelength))
        absorption_upper = np.zeros(len(wavelength))

        # Plot emission contributions
        for species_counter, identifier in enumerate(species):
            try:
                values = emission_df[tuple(identifier)].to_numpy()
                emission_upper = emission_lower + values
                color = color_list[species_counter]
                element_symbol = self._get_element_symbol(identifier[0])
                key = f"{element_symbol} (emission)"

                self.line_renderers[key] = p.varea(
                    x=wavelength,
                    y1=emission_lower,
                    y2=emission_upper,
                    fill_color=color,
                    fill_alpha=0.8,
                )

                # Create checkbox for this element emission
                checkbox = pn.widgets.Checkbox(name=key, value=True)
                checkbox.param.watch(self._on_toggle_line, "value")
                self.element_checkboxes[key] = checkbox

                # Store info for legend (only once per element)
                if not any(info[0] == element_symbol for info in self._element_info_for_legend):
                    self._element_info_for_legend.append((element_symbol, color, identifier[0]))

                emission_lower = emission_upper
            except KeyError:
                info_msg = (
                    f"{atomic_number2element_symbol(identifier)}"
                    f" is not in the emission packets; skipping"
                )
                logger.info(info_msg)
                pass

        # Plot absorption contributions
        for species_counter, identifier in enumerate(species):
            try:
                values = absorption_df[tuple(identifier)].to_numpy()
                absorption_lower = absorption_upper - values
                color = color_list[species_counter]
                element_symbol = self._get_element_symbol(identifier[0])
                key = f"{element_symbol} (absorption)"

                self.line_renderers[key] = p.varea(
                    x=wavelength,
                    y1=absorption_lower,
                    y2=absorption_upper,
                    fill_color=color,
                    fill_alpha=0.8,
                )

                # Create checkbox for this element absorption
                checkbox = pn.widgets.Checkbox(name=key, value=True)
                checkbox.param.watch(self._on_toggle_line, "value")
                self.element_checkboxes[key] = checkbox

                absorption_upper = absorption_lower
            except KeyError:
                info_msg = (
                    f"{atomic_number2element_symbol(identifier)}"
                    f" is not in the absorption packets; skipping"
                )
                logger.info(info_msg)
                pass

    def _add_element_line(self, p, key, x, y, color, label, visible=True, dash="solid", alpha=0.8):
        """Add an element contribution line to the plot."""
        self.line_renderers[key] = p.line(x, y, color=color, line_width=1.5,
                                          line_dash=dash, alpha=alpha, visible=visible)
        checkbox = pn.widgets.Checkbox(name=label, value=visible)
        checkbox.param.watch(self._on_toggle_line, "value")
        self.element_checkboxes[key] = checkbox

    def _create_element_legend_widget(self):
        """
        Create a vertical colorbar-style legend widget for elements using Panel.

        Returns
        -------
        pn.Column
            A Panel column containing the colorbar legend
        """
        if not hasattr(self, '_element_info_for_legend') or not self._element_info_for_legend:
            return None

        legend_items = []

        for element_symbol, color, atomic_num in self._element_info_for_legend:
            color_box = pn.pane.HTML(
                f"<div style='width: 30px; height: 25px; background-color: {color}; "
                f"border: 1px solid #ccc;'></div>",
                width=30, height=25
            )
            label = pn.pane.HTML(
                f"<span style='font-size: 12px; margin-left: 5px;'>{element_symbol}</span>",
                width=40
            )
            legend_items.append(pn.Row(color_box, label, height=27))

        legend_column = pn.Column(
            *legend_items,
            pn.pane.HTML("<b style='font-size: 11px;'>Elements</b>"),
            sizing_mode='fixed',
            width=80,
        )

        # Reverse order so highest contribution is at top
        legend_column = pn.Column(
            pn.pane.HTML("<b style='font-size: 11px;'>Elements</b>"),
            *reversed(legend_items),
            sizing_mode='fixed',
            width=80,
        )

        return legend_column

    def _get_element_symbol(self, atomic_num):
        """Get element symbol from atomic number."""
        try:
            return atomic_number2element_symbol(atomic_num)
        except (KeyError, ValueError):
            return f"Z={atomic_num}"

    def _selection_callback(self, _attr, _old, new):
        """
        Bokeh selection callback for spectrum plot.

        This method handles selection events from the Bokeh plot and updates
        the species interactions table based on the selected wavelength range.

        Parameters
        ----------
        _attr : str
            Attribute name (unused)
        _old : list
            Previous selection indices (unused)
        new : list
            New selection indices
        """
        if new:
            indices = new
            if len(indices) > 0:
                selected_x = [self._wavelength_data.value[i] for i in indices]
                x_range = [min(selected_x), max(selected_x)]

                # Update selection overlay to show persistent selection
                self._selection_source.data = dict(
                    left=[x_range[0]],
                    right=[x_range[1]],
                    top=[self._y_range[1]],
                    bottom=[self._y_range[0]],
                )

                # Track the current selection
                self._current_wavelength_range = x_range

                # Get current filter mode from buttons
                filter_mode_index = list(self.FILTER_MODES_DESC).index(
                    self.filter_mode_buttons.value
                )
                self._update_species_interactions(
                    x_range, self.FILTER_MODES[filter_mode_index]
                )

    def _update_species_interactions(self, wavelength_range, filter_mode):
        """
        Update data in species_interactions_table.

        The parameters are exact same as that of :code:`get_species_interactions`.
        Besides, it also does selection of 1st row in this table to trigger
        update in last_line_counts_table.
        """
        # Update data in species_interactions_table
        self.species_interactions_table.df = self.get_species_interactions(
            wavelength_range, filter_mode
        )

        # Get index of 0th row in species_interactions_table
        if (
            not self.species_interactions_table.df.empty
            and self.species_interactions_table.df.index[0] != ""
        ):
            species0 = self.species_interactions_table.df.index[0]

            # Also update last_line_counts_table by triggering its event listener
            if self.species_interactions_table.get_selected_rows() == [0]:
                # Listener won't trigger if last row selected in
                # species_interactions_table was also 0th, so unselect the rows
                self.species_interactions_table.change_selection([])
            # Select 0th row in this table to trigger _update_last_line_counts
            self.species_interactions_table.change_selection([species0])
        else:
            # Clear selection if no valid data
            self.species_interactions_table.change_selection([])

    def _update_last_line_counts(self, species, filter_mode, group_mode):
        """
        Update data in last_line_counts_table and associated total_packets_label.

        The parameters are exact same as that of :code:`get_last_line_counts`.
        """
        # Update data in line counts table
        self.last_line_counts_table.df = self.get_last_line_counts(
            species, filter_mode, group_mode
        )

        # Update its corresponding total_packets_label
        if species:
            self.total_packets_label.update_and_resize(
                self.last_line_counts_table.df.iloc[:, 0].sum()
            )
        else:  # Line counts table will be empty
            self.total_packets_label.update_and_resize(0)

    def _on_toggle_line(self, event):
        """Toggle visibility of spectrum lines based on widget state."""
        toggle_map = {
            "Real Packets": "real_packets",
            "Virtual Packets": "virtual_packets",
            "SDEC Spectrum": "sdec_spectrum",
            "Blackbody Photosphere": "photosphere",
            "No Interaction": "no_interaction",
            "Electron Scatter Only": "electron_scatter",
            "Observed Spectrum": "observed",
        }

        widget_name = event.obj.name
        renderer_key = toggle_map.get(widget_name)

        if renderer_key and renderer_key in self.line_renderers:
            self.line_renderers[renderer_key].visible = event.new
        else:
            for key, checkbox in self.element_checkboxes.items():
                if checkbox is event.obj:
                    if key in self.line_renderers:
                        self.line_renderers[key].visible = event.new
                    break

    def _filter_mode_toggle_handler(self, event):
        """
        Event handler for toggle in filter_mode_buttons.

        This method has the expected signature of the callback function
        for Panel widgets.
        """
        # Use tracked wavelength range instead of trying to access plotly shapes
        if self._current_wavelength_range is not None:
            # Get index from the selected value
            filter_mode_index = list(self.FILTER_MODES_DESC).index(event.new)
            self._update_species_interactions(
                self._current_wavelength_range,
                self.FILTER_MODES[filter_mode_index],
            )

    def _species_intrctn_selection_handler(self, event, _panel_widget):
        """
        Event handler for selection in species_interactions_table.

        This method has the expected signature of the function passed to
        :code:`handler` argument of :code:`on` method of PanelTableWidget.
        """
        # Don't execute function if no row was selected
        if not event["new"]:
            return

        # Get species from the selected row in species_interactions_table
        species_selected = self.species_interactions_table.df.index[
            event["new"][0]
        ]
        if species_selected == "":  # when species_interactions_table is empty
            species_selected = None

        # Get indices from the selected values
        filter_mode_index = list(self.FILTER_MODES_DESC).index(
            self.filter_mode_buttons.value
        )
        group_mode_index = list(self.GROUP_MODES_DESC).index(
            self.group_mode_dropdown.value
        )

        self._update_last_line_counts(
            species_selected,
            self.FILTER_MODES[filter_mode_index],
            self.GROUP_MODES[group_mode_index],
        )

    def _group_mode_dropdown_handler(self, event):
        """
        Event handler for selection in group_mode_dropdown.

        This method has the expected signature of the callback function
        for Panel widgets.
        """
        try:
            selected_row_idx = (
                self.species_interactions_table.get_selected_rows()[0]
            )
            species_selected = self.species_interactions_table.df.index[
                selected_row_idx
            ]
        except IndexError:  # No row is selected in species_interactions_table
            return

        # Get indices from the selected values
        filter_mode_index = list(self.FILTER_MODES_DESC).index(
            self.filter_mode_buttons.value
        )
        group_mode_index = list(self.GROUP_MODES_DESC).index(event.new)

        self._update_last_line_counts(
            species_selected,
            self.FILTER_MODES[filter_mode_index],
            self.GROUP_MODES[group_mode_index],
        )

    @staticmethod
    def ui_control_description(text):
        """Get description label of a UI control with increased font size."""
        return pn.pane.HTML(f"<span style='font-size: 1.15em;'>{text}:</span>")

    def display(self):
        """
        Display the fully-functional line info widget.

        It puts together all component widgets nicely together and enables
        interaction between all the components.

        Returns
        -------
        panel.Column
            Line info widget containing all component widgets
        """
        if not Environment.allows_widget_display():
            print("Please use a notebook to display the widget")
        else:
            # Panel tables handle their own sizing
            self.total_packets_label.update_and_resize(0)

            self.filter_mode_buttons.param.watch(
                self._filter_mode_toggle_handler, "value"
            )
            self.species_interactions_table.on(
                "selection_changed", self._species_intrctn_selection_handler
            )
            self.group_mode_dropdown.param.watch(
                self._group_mode_dropdown_handler, "value"
            )

            selection_box_symbol = (
                "<span style='display: inline-block; "
                f"background-color: {self.COLORS['selection_area']}; "
                f"border: 1px solid {self.COLORS['selection_border']}; "
                "width: 0.8em; height: 1.2em; vertical-align: middle;'></span>"
            )

            visibility_checkboxes = [
                self.checkbox_real_packets,
                self.checkbox_virtual_packets,
            ]
            if self.show_sdec:
                visibility_checkboxes.extend([
                    self.checkbox_sdec_spectrum,
                    self.checkbox_photosphere,
                    self.checkbox_no_interaction,
                    self.checkbox_electron_scatter,
                ])
            if self.observed_spectrum is not None:
                visibility_checkboxes.append(self.checkbox_observed_spectrum)

            row_size = 5
            visibility_rows = [
                pn.Row(*visibility_checkboxes[i:i + row_size])
                for i in range(0, len(visibility_checkboxes), row_size)
            ]

            visibility_panel = pn.Card(
                *visibility_rows,
                title="Line Visibility",
                collapsed=False,
            )

            element_panel = None
            if self.show_sdec and self.element_checkboxes:
                element_checkbox_list = list(self.element_checkboxes.values())
                element_rows = [
                    pn.Row(*element_checkbox_list[i:i + row_size])
                    for i in range(0, len(element_checkbox_list), row_size)
                ]
                element_panel = pn.Card(
                    *element_rows,
                    title="Species Contributions",
                    collapsed=True,
                )

            # Create Panel description components
            filter_description = pn.pane.HTML(
                f"<span style='font-size: 1.15em;'>Filter selected wavelength range "
                f"( {selection_box_symbol} ) by:</span>"
            )

            group_description = pn.pane.HTML(
                "<span style='font-size: 1.15em;'>Group packet counts by:</span>"
            )

            table_container_left = pn.Column(
                filter_description,
                self.filter_mode_buttons,
                self.species_interactions_table.table,
                margin=(0, 15),
            )

            table_container_right = pn.Column(
                group_description,
                self.group_mode_dropdown,
                self.last_line_counts_table.table,
                self.total_packets_label.widget,
                margin=(0, 15),
            )

            tables_row = pn.Row(
                table_container_left,
                table_container_right,
                sizing_mode="stretch_width",
            )

            element_legend = None
            if self.show_sdec and hasattr(self, '_element_info_for_legend'):
                element_legend = self._create_element_legend_widget()

            # Arrange figure with element legend on the right
            if element_legend is not None:
                figure_with_legend = pn.Row(
                    self.figure_widget,
                    element_legend,
                    sizing_mode="stretch_width"
                )
            else:
                figure_with_legend = self.figure_widget

            panels = [visibility_panel]
            if element_panel is not None:
                panels.append(element_panel)
            panels.extend([figure_with_legend, tables_row])

            widget = pn.Column(
                *panels, sizing_mode="stretch_width"
            )

            return widget
