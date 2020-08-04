from astropy import units as u
import pandas as pd
import qgrid
import plotly.graph_objects as go

from tardis.analysis import LastLineInteraction
from tardis.util.base import species_tuple_to_string, species_string_to_tuple


class LineInfoWidget:
    def __init__(
        self,
        lines_data,
        line_interaction_analysis,
        spectrum_wavelength,
        spectrum_luminosity_density_lambda,
    ):
        self.lines_data = lines_data
        self.line_interaction_analysis = line_interaction_analysis

        # Widgets
        self.species_abundances_table = self.create_table_widget(
            self.get_species_abundances(None), [35, 65]
        )
        self.line_counts_table = self.create_table_widget(
            self.get_last_line_counts(None), [75, 25]
        )
        self.figure_widget = self.plot_spectrum(
            spectrum_wavelength, spectrum_luminosity_density_lambda
        )

    @classmethod
    def from_simulation(cls, sim):
        return cls(
            lines_data=sim.plasma.lines.reset_index().set_index("line_id"),
            line_interaction_analysis=LastLineInteraction.from_model(sim),
            spectrum_wavelength=sim.runner.spectrum.wavelength,
            spectrum_luminosity_density_lambda=sim.runner.spectrum.luminosity_density_lambda,
        )

    def get_species_abundances(
        self, wavelength_range, packet_filter_mode="packet_nu"
    ):
        if wavelength_range:
            self.line_interaction_analysis.packet_filter_mode = (
                packet_filter_mode
            )
            self.line_interaction_analysis.wavelength_start = (
                wavelength_range[0] * u.angstrom
            )
            self.line_interaction_analysis.wavelength_end = (
                wavelength_range[1] * u.angstrom
            )

            selected_species_group = self.line_interaction_analysis.last_line_in.groupby(
                ["atomic_number", "ion_number"]
            )
            selected_species_symbols = [
                species_tuple_to_string(item)
                for item in selected_species_group.groups.keys()
            ]

            selected_species_abundances = (
                selected_species_group.size()
                / self.line_interaction_analysis.last_line_in.shape[0]
            )

        else:
            # Qgrid can't create index col of table widget with empty dataframes
            # so create one row with null strings to emulate empty table widget
            selected_species_symbols = [""]
            selected_species_abundances = pd.Series([""])

        selected_species_abundances.index = pd.Index(
            selected_species_symbols, name="Species"
        )
        selected_species_abundances.name = "Fractional Abundance"
        return selected_species_abundances.sort_values(
            ascending=False
        ).to_frame()

    def get_last_line_counts(self, selected_species):
        if selected_species:
            selected_species_tuple = species_string_to_tuple(selected_species)

            current_last_line_in = self.line_interaction_analysis.last_line_in.xs(
                key=(selected_species_tuple[0], selected_species_tuple[1]),
                level=["atomic_number", "ion_number"],
                drop_level=False,
            ).reset_index()

            current_last_line_out = self.line_interaction_analysis.last_line_out.xs(
                key=(selected_species_tuple[0], selected_species_tuple[1]),
                level=["atomic_number", "ion_number"],
                drop_level=False,
            ).reset_index()

            current_last_line_in["line_id_out"] = current_last_line_out.line_id

            last_line_in_string = []
            last_line_count = []
            grouped_line_interactions = current_last_line_in.groupby(
                ["line_id", "line_id_out"]
            )

            for line_id, count in grouped_line_interactions.size().iteritems():
                current_line_in = self.lines_data.loc[line_id[0]]
                current_line_out = self.lines_data.loc[line_id[1]]
                last_line_in_string.append(
                    f"exc. {int(current_line_in.level_number_lower)}-"
                    f"{int(current_line_in.level_number_upper)} "
                    f"({current_line_in.wavelength:.2f} A) "
                    f"de-exc. {int(current_line_out.level_number_upper)}-"
                    f"{int(current_line_out.level_number_lower)} "
                    f"({current_line_out.wavelength:.2f} A)"
                )
                last_line_count.append(count)

        else:
            last_line_count = [""]
            last_line_in_string = [""]

        line_counts = pd.Series(last_line_count)
        line_counts.name = "No. of packets"
        line_counts.index = pd.Index(
            last_line_in_string, name="Last Line Interaction"
        )
        return line_counts.sort_values(ascending=False).to_frame()

    @staticmethod
    def create_table_widget(data, col_widths):
        grid_options = {
            "sortable": False,
            "filterable": False,
            "editable": False,
            "minVisibleRows": 2,
            "maxVisibleRows": 10,
        }
        column_options = {
            "minWidth": None,
        }

        # Preparing dictionary that defines column widths
        cols_with_index = [data.index.name] + data.columns.to_list()
        column_widths_definitions = {
            col_name: {"width": col_width}
            for col_name, col_width in zip(cols_with_index, col_widths)
        }

        return qgrid.show_grid(
            data,
            grid_options=grid_options,
            column_options=column_options,
            column_definitions=column_widths_definitions,
        )

    @staticmethod
    def axis_label_in_latex(label_text, unit):
        unit_in_latex = unit.to_string("latex_inline").strip("$")
        return f"$\\text{{{label_text}}}\\,[{unit_in_latex}]$"

    def plot_spectrum(self, wavelength, luminosity_density_lambda):
        initial_zoomed_range = [
            3000,
            9000,
        ]  # TODO: select from wavelength logically

        return go.FigureWidget(
            [
                go.Scatter(
                    x=wavelength, y=luminosity_density_lambda, showlegend=False
                ),
                # Hide a one-point scatter trace, to bring boxselect in modebar
                go.Scatter(
                    x=wavelength[0],
                    y=luminosity_density_lambda[0],
                    mode="markers",
                    marker=dict(opacity=0),
                    showlegend=False,
                ),
            ],
            layout=go.Layout(
                title="Spectrum of real packets",
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
            ),
        )

