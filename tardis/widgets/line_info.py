from astropy import units as u
import numpy as np
import pandas as pd
import qgrid
from plotly import graph_objects as go
from plotly.callbacks import BoxSelector
import ipywidgets as ipw

from tardis.analysis import LastLineInteraction
from tardis.util.base import species_tuple_to_string, species_string_to_tuple
from tardis.widgets.util import create_table_widget, TableSummaryLabel


class LineInfoWidget:
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
        self.lines_data = lines_data
        self.line_interaction_analysis = line_interaction_analysis

        # Widgets ---------------------------------------------
        max_rows_option = {"maxVisibleRows": 9}
        self.species_abundances_table = create_table_widget(
            data=self.get_species_abundances(None),
            col_widths=[35, 65],
            table_options=max_rows_option,
        )

        line_counts_col_widths = [75, 25]
        self.line_counts_table = create_table_widget(
            data=self.get_last_line_counts(None),
            col_widths=line_counts_col_widths,
            table_options=max_rows_option,
        )
        self.total_packets_label = TableSummaryLabel(
            target_table=self.line_counts_table,
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
        return cls(
            lines_data=sim.plasma.lines.reset_index().set_index("line_id"),
            line_interaction_analysis={
                filter_mode: LastLineInteraction.from_model(sim, filter_mode)
                for filter_mode in cls.FILTER_MODES
            },
            spectrum_wavelength=sim.runner.spectrum.wavelength,
            spectrum_luminosity_density_lambda=sim.runner.spectrum.luminosity_density_lambda,
            virt_spectrum_wavelength=sim.runner.spectrum_virtual.wavelength,
            virt_spectrum_luminosity_density_lambda=sim.runner.spectrum_virtual.luminosity_density_lambda,
        )

    def get_species_abundances(
        self, wavelength_range, filter_mode=FILTER_MODES[0]
    ):
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

                # Normalize each group's count to find fractional abundances
                selected_species_abundances = (
                    selected_species_group.size()
                    / self.line_interaction_analysis[
                        filter_mode
                    ].last_line_in.shape[0]
                )

            else:
                selected_species_symbols = [""]
                selected_species_abundances = pd.Series([""])

        else:  # wavelength_range is None or ""
            selected_species_symbols = [""]
            selected_species_abundances = pd.Series([""])

        selected_species_abundances.index = pd.Index(
            selected_species_symbols, name="Species"
        )
        selected_species_abundances.name = "Fraction of packets interacting"
        return selected_species_abundances.sort_values(
            ascending=False
        ).to_frame()

    def get_last_line_counts(
        self,
        selected_species,
        filter_mode=FILTER_MODES[0],
        group_mode=GROUP_MODES[0],
    ):
        if selected_species:
            selected_species_tuple = species_string_to_tuple(selected_species)

            # Get selected species' rows from last_line_in dataframe
            current_last_lines_in = (
                self.line_interaction_analysis[filter_mode]
                .last_line_in.xs(
                    key=(selected_species_tuple[0], selected_species_tuple[1]),
                    level=["atomic_number", "ion_number"],
                    drop_level=False,
                )
                .reset_index()
            )

            # Get selected species' rows from last_line_out dataframe
            current_last_lines_out = (
                self.line_interaction_analysis[filter_mode]
                .last_line_out.xs(
                    key=(selected_species_tuple[0], selected_species_tuple[1]),
                    level=["atomic_number", "ion_number"],
                    drop_level=False,
                )
                .reset_index()
            )

            last_line_interaction_string = []
            last_line_count = []

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
                ) in grouped_line_interactions.size().iteritems():
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
                    last_line_count.append(count)

            elif group_mode == "exc":
                grouped_line_interactions = current_last_lines_in.groupby(
                    "line_id"
                )

                # Iterate over each group's key and size and append them to list
                for (
                    line_id,
                    count,
                ) in grouped_line_interactions.size().iteritems():
                    current_line_in = self.lines_data.loc[line_id]
                    last_line_interaction_string.append(
                        f"exc. {int(current_line_in.level_number_lower):02d}-"
                        f"{int(current_line_in.level_number_upper):02d} "
                        f"({current_line_in.wavelength:.2f} A)"
                    )
                    last_line_count.append(count)

            elif group_mode == "de-exc":
                grouped_line_interactions = current_last_lines_out.groupby(
                    "line_id"
                )

                # Iterate over each group's key and size and append them to list
                for (
                    line_id,
                    count,
                ) in grouped_line_interactions.size().iteritems():
                    current_line_out = self.lines_data.loc[line_id]
                    last_line_interaction_string.append(
                        f"de-exc. {int(current_line_out.level_number_upper):02d}-"
                        f"{int(current_line_out.level_number_lower):02d} "
                        f"({current_line_out.wavelength:.2f} A)"
                    )
                    last_line_count.append(count)

            else:
                raise ValueError(
                    "Invalid value passed to group_mode argument. "
                    f"Allowed values are {self.GROUP_MODES}"
                )

        else:
            last_line_count = [""]
            last_line_interaction_string = [""]

        line_counts = pd.Series(last_line_count)
        line_counts.name = "No. of packets"
        line_counts.index = pd.Index(
            last_line_interaction_string, name="Last Line Interaction"
        )
        return line_counts.sort_values(ascending=False).to_frame()

    @staticmethod
    def axis_label_in_latex(label_text, unit):
        # TODO: If erg and s-1 present stick them together as erg and s
        unit_in_latex = unit.to_string("latex_inline").strip("$")
        return f"$\\text{{{label_text}}}\\,[{unit_in_latex}]$"

    @staticmethod
    def get_middle_half_edges(arr):
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
        initial_zoomed_range = self.get_middle_half_edges(wavelength.value)

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
                    x=wavelength[0],
                    y=luminosity_density_lambda[0],
                    mode="markers",
                    marker=dict(opacity=0),
                    showlegend=False,
                ),
            ],
            layout=go.Layout(
                title="Spectrum",
                xaxis=dict(
                    title=self.axis_label_in_latex(
                        "Wavelength", wavelength.unit
                    ),
                    exponentformat="none",
                    rangeslider=dict(visible=True),
                    range=initial_zoomed_range,
                ),
                yaxis=dict(
                    title=self.axis_label_in_latex(
                        "Luminosity",
                        luminosity_density_lambda.to("erg/(s AA)").unit,
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

    def update_species_abundances(self, wavelength_range, filter_mode):
        # Update data in species abundance table
        self.species_abundances_table.df = self.get_species_abundances(
            wavelength_range, filter_mode
        )

        # Get index of 0th row in species abundance table
        species0 = self.species_abundances_table.df.index[0]

        # Also update line counts table by triggering its event listener
        # Listener won't trigger if last row selected in species abundance table was also 0th
        if self.species_abundances_table.get_selected_rows() == [0]:
            self.species_abundances_table.change_selection([])  # Unselect rows
        # Select 0th row in count table which will trigger update_last_line_counts
        self.species_abundances_table.change_selection([species0])

    def add_selection_box(self, selector):
        self.figure_widget.layout.shapes = [
            dict(
                type="rect",
                xref="x",
                yref="y",
                x0=selector.xrange[0],
                y0=selector.yrange[0],
                x1=selector.xrange[1],
                y1=selector.yrange[1],
                line=dict(color=self.COLORS["selection_border"], width=1,),
                fillcolor=self.COLORS["selection_area"],
                opacity=0.5,
            )
        ]

    def update_last_line_counts(self, species, filter_mode, group_mode):
        # Update data in line counts table
        self.line_counts_table.df = self.get_last_line_counts(
            species, filter_mode, group_mode
        )

        # Update its corresponding total packets label
        if species:
            self.total_packets_label.update_and_resize(
                self.line_counts_table.df.iloc[:, 0].sum()
            )
        else:  # Line counts table will be empty
            self.total_packets_label.update_and_resize(0)

    def spectrum_selection_handler(self, trace, points, selector):
        if isinstance(selector, BoxSelector):
            self.add_selection_box(selector)
            self.update_species_abundances(
                selector.xrange,
                self.FILTER_MODES[self.filter_mode_buttons.index],
            )

    def filter_mode_toggle_handler(self, change):
        try:
            wavelength_range = [
                self.figure_widget.layout.shapes[0][x] for x in ("x0", "x1")
            ]
        except IndexError:  # No selection is made on figure widget
            return

        self.update_species_abundances(
            wavelength_range, self.FILTER_MODES[self.filter_mode_buttons.index],
        )

    def species_abund_selection_handler(self, event, qgrid_widget):
        # Don't execute function if no row was selected implicitly (by api)
        if event["new"] == [] and event["source"] == "api":
            return

        # Get species from selected row in species abundance table
        species_selected = self.species_abundances_table.df.index[
            event["new"][0]
        ]

        self.update_last_line_counts(
            species_selected,
            self.FILTER_MODES[self.filter_mode_buttons.index],
            self.GROUP_MODES[self.group_mode_dropdown.index],
        )

    def group_mode_dropdown_handler(self, change):
        try:
            selected_row_idx = self.species_abundances_table.get_selected_rows()[
                0
            ]
            species_selected = self.species_abundances_table.df.index[
                selected_row_idx
            ]
        except IndexError:  # No row is selected in species abundances table
            return

        self.update_last_line_counts(
            species_selected,
            self.FILTER_MODES[self.filter_mode_buttons.index],
            self.GROUP_MODES[self.group_mode_dropdown.index],
        )

    @staticmethod
    def ui_control_description(text):
        return ipw.HTML(f"<span style='font-size: 1.15em;'>{text}:</span>")

    def display(self):
        # Set widths of widgets
        self.species_abundances_table.layout.width = "350px"
        self.line_counts_table.layout.width = "450px"
        self.total_packets_label.update_and_resize(0)
        self.group_mode_dropdown.layout.width = "auto"

        # Attach event listeners to widgets
        spectrum_trace = self.figure_widget.data[0]
        spectrum_trace.on_selection(self.spectrum_selection_handler)
        self.filter_mode_buttons.observe(
            self.filter_mode_toggle_handler, names="index"
        )
        self.species_abundances_table.on(
            "selection_changed", self.species_abund_selection_handler
        )
        self.group_mode_dropdown.observe(
            self.group_mode_dropdown_handler, names="index"
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
                self.species_abundances_table,
            ],
            layout=dict(margin="0px 15px"),
        )

        table_container_right = ipw.VBox(
            [
                self.ui_control_description("Group packet counts by"),
                self.group_mode_dropdown,
                self.line_counts_table,
                self.total_packets_label.widget,
            ],
            layout=dict(margin="0px 15px"),
        )

        return ipw.VBox(
            [
                self.figure_widget,
                ipw.Box(
                    [table_container_left, table_container_right,],
                    layout=dict(
                        display="flex",
                        align_items="flex-start",
                        justify_content="center",
                    ),
                ),
            ]
        )
