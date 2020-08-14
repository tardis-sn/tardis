from astropy import units as u
import numpy as np
import pandas as pd
import qgrid
import plotly.graph_objects as go
import ipywidgets as ipw

from tardis.analysis import LastLineInteraction
from tardis.util.base import species_tuple_to_string, species_string_to_tuple
from tardis.widgets.util import create_table_widget, TableSummaryLabel


class LineInfoWidget:
    filter_modes = ["packet_out_nu", "packet_in_nu"]
    filter_modes_description = ["Emitted Wavelength", "Absorbed Wavelength"]
    group_modes = ["both", "exc", "de-exc"]
    group_modes_descripton = [
        "Both excitation line (absorption) and de-excitation line (emission)",
        "Only excitation line (absorption)",
        "Only de-excitation line (emission)",
    ]

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
            options=self.filter_modes_description, index=0
        )

        self.group_mode_dropdown = ipw.Dropdown(
            options=self.group_modes_descripton,
            index=0,
            layout=dict(
                width="auto"
            ),  # So that it can take all width available
        )

    @classmethod
    def from_simulation(cls, sim):
        return cls(
            lines_data=sim.plasma.lines.reset_index().set_index("line_id"),
            line_interaction_analysis={
                filter_mode: LastLineInteraction.from_model(sim, filter_mode)
                for filter_mode in cls.filter_modes
            },
            spectrum_wavelength=sim.runner.spectrum.wavelength,
            spectrum_luminosity_density_lambda=sim.runner.spectrum.luminosity_density_lambda,
            virt_spectrum_wavelength=sim.runner.spectrum_virtual.wavelength,
            virt_spectrum_luminosity_density_lambda=sim.runner.spectrum_virtual.luminosity_density_lambda,
        )

    def get_species_abundances(
        self, wavelength_range, filter_mode=filter_modes[0]
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

        else:
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
        filter_mode=filter_modes[0],
        group_mode=group_modes[0],
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
                    f"Allowed values are {self.group_modes}"
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

        figure_widget = go.FigureWidget(
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

        return figure_widget
