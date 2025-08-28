"""Class to create and display Line Info Widget."""

import numpy as np
import pandas as pd
import panel as pn
from astropy import units as u

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from tardis.analysis import LastLineInteraction
from tardis.util.base import (
    species_string_to_tuple,
    species_tuple_to_string,
)

from tardis.visualization.widgets.util import (
    TableSummaryLabel,
    create_table_widget,
)
from tardis.util.environment import Environment


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
        """
        self.lines_data = lines_data
        self.line_interaction_analysis = line_interaction_analysis

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
            options=list(self.FILTER_MODES_DESC), value=self.FILTER_MODES_DESC[0]
        )

        self.group_mode_dropdown = pn.widgets.Select(
            options=list(self.GROUP_MODES_DESC), value=self.GROUP_MODES_DESC[0]
        )

        self._current_wavelength_range = None  # Track current selection

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of LineInfoWidget from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        LineInfoWidget object
        """
        spectrum_solver = sim.spectrum_solver
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
            ascending=False
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
                current_last_lines_in[
                    "line_id_out"
                ] = current_last_lines_out.line_id
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
        return last_line_counts.sort_values(ascending=False).to_frame()

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
        arr = np.sort(arr)
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
        Produce a plotly figure widget by plotting the spectrum of model.

        Parameters
        ----------
        wavelength : astropy.Quantity
            Wavelength values of a real spectrum, having unit of Angstrom
        luminosity_density_lambda : astropy.Quantity
            Luminosity density lambda values of a real spectrum, having unit
            of (erg/s)/Angstrom
        virt_wavelength : astropy.Quantity
            Wavelength values of a virtual spectrum, having unit of Angstrom
        virt_luminosity_density_lambda : astropy.Quantity
            Luminosity density lambda values of a virtual spectrum, having unit
            of (erg/s)/Angstrom

        Returns
        -------
        plotly.graph_objects.FigureWidget
        """
        # Create Bokeh figure with box select
        p = figure(width=800, height=400, title='Spectrum', tools='box_select,reset')
        
        # Add proper axis labels with units (Bokeh compatible)
        p.xaxis.axis_label = f"Wavelength [{wavelength.unit}]"
        p.yaxis.axis_label = f"Luminosity [{luminosity_density_lambda.unit}]"
        
        # Add line plots
        p.line(wavelength.value, luminosity_density_lambda.value, legend_label='Real packets', color='blue')
        p.line(virt_wavelength.value, virt_luminosity_density_lambda.value, legend_label='Virtual packets', color='red')
        
        # Create invisible scatter for selection (needed for box select to work)
        source = ColumnDataSource(dict(x=wavelength.value, y=luminosity_density_lambda.value))
        p.scatter('x', 'y', source=source, alpha=0, size=1)
        
        # Create selection overlay source (initially empty)
        selection_source = ColumnDataSource(dict(left=[], right=[], top=[], bottom=[]))
        p.quad(left='left', right='right', top='top', bottom='bottom',
               source=selection_source, alpha=0.3, color='lightblue')
        
        # Store references for callback
        self._bokeh_plot = p
        self._selection_source = selection_source
        self._y_range = [luminosity_density_lambda.value.min(), luminosity_density_lambda.value.max()]
        self._wavelength_data = wavelength  # Store wavelength data for callback access

        # Connect selection callback
        source.selected.on_change('indices', self._selection_callback)
        
        return pn.pane.Bokeh(p)

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
                    bottom=[self._y_range[0]]
                )

                # Track the current selection
                self._current_wavelength_range = x_range

                # Get current filter mode from buttons
                filter_mode_index = list(self.FILTER_MODES_DESC).index(self.filter_mode_buttons.value)
                self._update_species_interactions(x_range, self.FILTER_MODES[filter_mode_index])

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
        if not self.species_interactions_table.df.empty and self.species_interactions_table.df.index[0] != "":
            
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
        filter_mode_index = list(self.FILTER_MODES_DESC).index(self.filter_mode_buttons.value)
        group_mode_index = list(self.GROUP_MODES_DESC).index(self.group_mode_dropdown.value)

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
        filter_mode_index = list(self.FILTER_MODES_DESC).index(self.filter_mode_buttons.value)
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
                margin=(0, 15)
            )

            table_container_right = pn.Column(
                group_description,
                self.group_mode_dropdown,
                self.last_line_counts_table.table,
                self.total_packets_label.widget,
                margin=(0, 15)
            )

            tables_row = pn.Row(
                table_container_left,
                table_container_right,
                sizing_mode='stretch_width'
            )

            widget = pn.Column(
                self.figure_widget,
                tables_row,
                sizing_mode='stretch_width'
            )
            
            return widget

