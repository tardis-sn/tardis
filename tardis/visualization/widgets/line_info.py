"""Class to create and display Line Info Widget."""
import re

from astropy import units as u
import numpy as np
import pandas as pd
import qgrid
from plotly import graph_objects as go
from plotly.callbacks import BoxSelector
import ipywidgets as ipw

from tardis.analysis import LastLineInteraction
from tardis.util.base import (
    species_tuple_to_string,
    species_string_to_tuple,
    is_notebook,
)
from tardis.visualization.widgets.util import (
    create_table_widget,
    TableSummaryLabel,
)
from tardis.visualization import plot_util as pu


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
            col_widths=[35, 65],
            table_options=max_rows_option,
        )

        line_counts_col_widths = [75, 25]
        self.last_line_counts_table = create_table_widget(
            data=self.get_last_line_counts(None),
            col_widths=line_counts_col_widths,
            table_options=max_rows_option,
        )
        self.total_packets_label = TableSummaryLabel(
            target_table=self.last_line_counts_table,
            table_col_widths=line_counts_col_widths,
            label_key="Total Packets",
            label_value=0,
        )

        self.figure_widget = self.plot_spectrum(
            spectrum_wavelength,
            spectrum_luminosity_density_lambda,
            virt_spectrum_wavelength,
            virt_spectrum_luminosity_density_lambda,
        )

        self.filter_mode_buttons = ipw.ToggleButtons(
            options=self.FILTER_MODES_DESC, index=0
        )

        self.group_mode_dropdown = ipw.Dropdown(
            options=self.GROUP_MODES_DESC, index=0
        )

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
        transport_state = sim.transport.transport_state
        return cls(
            lines_data=sim.plasma.lines.reset_index().set_index("line_id"),
            line_interaction_analysis={
                filter_mode: LastLineInteraction.from_simulation(
                    sim, filter_mode
                )
                for filter_mode in cls.FILTER_MODES
            },
            spectrum_wavelength=transport_state.spectrum.wavelength,
            spectrum_luminosity_density_lambda=transport_state.spectrum.luminosity_density_lambda.to(
                "erg/(s AA)"
            ),
            virt_spectrum_wavelength=transport_state.spectrum_virtual.wavelength,
            virt_spectrum_luminosity_density_lambda=transport_state.spectrum_virtual.luminosity_density_lambda.to(
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
                # qgrid cannot show empty dataframe properly,
                # so create one row with empty strings
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
                    == False
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
            # qgrid cannot show empty dataframe properly,
            # so create one row with empty strings
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
        # Initially zoomed range in rangeslider should be middle half of spectrum
        initial_zoomed_range = self.get_middle_half_edges(wavelength.value)

        # The scatter point should be a middle point in spectrum otherwise
        # the extra padding around it will be oddly visible when near the edge
        scatter_point_idx = pu.get_mid_point_idx(wavelength.value)

        return go.FigureWidget(
            [
                go.Scatter(
                    x=wavelength,
                    y=luminosity_density_lambda,
                    name="Real packets",
                ),
                go.Scatter(
                    x=virt_wavelength,
                    y=virt_luminosity_density_lambda,
                    name="Virtual packets",
                ),
                # Hide a one point scatter trace, to bring boxselect in modebar
                go.Scatter(
                    x=wavelength[scatter_point_idx],
                    y=luminosity_density_lambda[scatter_point_idx],
                    mode="markers",
                    marker=dict(opacity=0),
                    showlegend=False,
                ),
            ],
            layout=go.Layout(
                title="Spectrum",
                xaxis=dict(
                    title=pu.axis_label_in_latex("Wavelength", wavelength.unit),
                    exponentformat="none",
                    rangeslider=dict(visible=True),
                    range=initial_zoomed_range,
                ),
                yaxis=dict(
                    title=pu.axis_label_in_latex(
                        "Luminosity",
                        luminosity_density_lambda.unit,
                    ),
                    exponentformat="e",
                    fixedrange=False,
                ),
                dragmode="select",
                selectdirection="h",
                height=400,
                margin=dict(t=50, b=60),
            ),
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
        species0 = self.species_interactions_table.df.index[0]

        # Also update last_line_counts_table by triggering its event listener
        if self.species_interactions_table.get_selected_rows() == [0]:
            # Listener won't trigger if last row selected in
            # species_interactions_table was also 0th, so unselect the rows
            self.species_interactions_table.change_selection([])
        # Select 0th row in this table to trigger _update_last_line_counts
        self.species_interactions_table.change_selection([species0])

    def _add_selection_box(self, selector):
        """
        Draw a shape on plotly figure widget to represent the selection.

        Parameters
        ----------
        selector : plotly.callbacks.BoxSelector
            The object containing data about current selection made on plot
            (x-axis and y-axis range of selection box)
        """
        self.figure_widget.layout.shapes = [
            dict(
                type="rect",
                xref="x",
                yref="y",
                x0=selector.xrange[0],
                y0=selector.yrange[0],
                x1=selector.xrange[1],
                y1=selector.yrange[1],
                line=dict(
                    color=self.COLORS["selection_border"],
                    width=1,
                ),
                fillcolor=self.COLORS["selection_area"],
                opacity=0.5,
            )
        ]

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

    def _spectrum_selection_handler(self, trace, points, selector):
        """
        Event handler for selection of spectrum in plotly figure widget.

        This method has the expected signature of the callback function passed
        to :code:`on_selection` method of a plotly trace as explained in
        `their docs <https://plotly.com/python-api-reference/generated/plotly.html#plotly.basedatatypes.BaseTraceType.on_selection>`_.
        """
        if isinstance(selector, BoxSelector):
            self._add_selection_box(selector)
            self._update_species_interactions(
                selector.xrange,
                self.FILTER_MODES[self.filter_mode_buttons.index],
            )

    def _filter_mode_toggle_handler(self, change):
        """
        Event handler for toggle in filter_mode_buttons.

        This method has the expected signature of the callback function
        passed to :code:`observe` method of ipywidgets as explained in
        `their docs <https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html#Signatures>`_.
        """
        try:
            wavelength_range = [
                self.figure_widget.layout.shapes[0][x] for x in ("x0", "x1")
            ]
        except IndexError:  # No selection is made on figure widget
            return

        self._update_species_interactions(
            wavelength_range,
            self.FILTER_MODES[self.filter_mode_buttons.index],
        )

    def _species_intrctn_selection_handler(self, event, qgrid_widget):
        """
        Event handler for selection in species_interactions_table.

        This method has the expected signature of the function passed to
        :code:`handler` argument of :code:`on_selection` method of qgrid.QgridWidget
        as explained in `their docs <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.on>`_.
        """
        # Don't execute function if no row was selected implicitly (by api)
        if event["new"] == [] and event["source"] == "api":
            return

        # Get species from the selected row in species_interactions_table
        species_selected = self.species_interactions_table.df.index[
            event["new"][0]
        ]
        if species_selected == "":  # when species_interactions_table is empty
            species_selected = None

        self._update_last_line_counts(
            species_selected,
            self.FILTER_MODES[self.filter_mode_buttons.index],
            self.GROUP_MODES[self.group_mode_dropdown.index],
        )

    def _group_mode_dropdown_handler(self, change):
        """
        Event handler for selection in group_mode_dropdown.

        This method has the expected signature of the callback function
        passed to :code:`observe` method of ipywidgets as explained in
        `their docs <https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html#Signatures>`_.
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

        self._update_last_line_counts(
            species_selected,
            self.FILTER_MODES[self.filter_mode_buttons.index],
            self.GROUP_MODES[self.group_mode_dropdown.index],
        )

    @staticmethod
    def ui_control_description(text):
        """Get description label of a UI control with increased font size."""
        return ipw.HTML(f"<span style='font-size: 1.15em;'>{text}:</span>")

    def display(self):
        """
        Display the fully-functional line info widget.

        It puts together all component widgets nicely together and enables
        interaction between all the components.

        Returns
        -------
        ipywidgets.Box
            Line info widget containing all component widgets
        """
        if not is_notebook():
            print("Please use a notebook to display the widget")
        else:
            # Set widths of widgets
            self.species_interactions_table.layout.width = "350px"
            self.last_line_counts_table.layout.width = "450px"
            self.total_packets_label.update_and_resize(0)
            self.group_mode_dropdown.layout.width = "auto"

            # Attach event listeners to widgets
            spectrum_trace = self.figure_widget.data[0]
            spectrum_trace.on_selection(self._spectrum_selection_handler)
            self.filter_mode_buttons.observe(
                self._filter_mode_toggle_handler, names="index"
            )
            self.species_interactions_table.on(
                "selection_changed", self._species_intrctn_selection_handler
            )
            self.group_mode_dropdown.observe(
                self._group_mode_dropdown_handler, names="index"
            )

            selection_box_symbol = (
                "<span style='display: inline-block; "
                f"background-color: {self.COLORS['selection_area']}; "
                f"border: 1px solid {self.COLORS['selection_border']}; "
                "width: 0.8em; height: 1.2em; vertical-align: middle;'></span>"
            )

            table_container_left = ipw.VBox(
                [
                    self.ui_control_description(
                        "Filter selected wavelength range "
                        f"( {selection_box_symbol} ) by"
                    ),
                    self.filter_mode_buttons,
                    self.species_interactions_table,
                ],
                layout=dict(margin="0px 15px"),
            )

            table_container_right = ipw.VBox(
                [
                    self.ui_control_description("Group packet counts by"),
                    self.group_mode_dropdown,
                    self.last_line_counts_table,
                    self.total_packets_label.widget,
                ],
                layout=dict(margin="0px 15px"),
            )

            return ipw.VBox(
                [
                    self.figure_widget,
                    ipw.Box(
                        [
                            table_container_left,
                            table_container_right,
                        ],
                        layout=dict(
                            display="flex",
                            align_items="flex-start",
                            justify_content="center",
                            height="420px",
                        ),
                    ),
                ]
            )
