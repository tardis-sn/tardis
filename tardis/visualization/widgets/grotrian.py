from tardis.analysis import LastLineInteraction
from tardis.util.base import int_to_roman
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
import ipywidgets as ipw


class GrotrianWidget:
    FILTER_MODES = ("packet_out_nu", "packet_in_nu")
    FILTER_MODES_DESC = ("Emitted Wavelength", "Absorbed Wavelength")

    @classmethod
    def from_simulation(cls, sim, **kwargs):
        atom_data = sim.plasma.atomic_data.atom_data
        level_energy_data = pd.Series(
            sim.plasma.atomic_data.levels.energy
            / (u.electronvolt.to("J") * 1e7),
            name="energy",
        )
        level_population_data = sim.plasma.level_number_density
        line_interaction_analysis = {
            filter_mode: LastLineInteraction.from_model(sim, filter_mode)
            for filter_mode in cls.FILTER_MODES
        }
        return cls(
            atom_data=atom_data,
            level_energy_data=level_energy_data,
            level_population_data=level_population_data,
            line_interaction_analysis=line_interaction_analysis,
            **kwargs,
        )

    def __init__(
        self,
        atom_data,
        level_energy_data,
        level_population_data,
        line_interaction_analysis,
        colorscale="plasma",
    ):
        # Set data members
        self.atom_data = atom_data
        self.level_energy_data = level_energy_data
        self.level_population_data = level_population_data
        self.line_interaction_analysis = line_interaction_analysis

        # Selector for max level threshold
        self.max_level_selector = ipw.BoundedIntText(
            value=10,
            min=1,
            max=40,
            step=1,
            description="Max Levels:",
        )

        # Energy difference threshold below which levels are merged
        self.level_diff_threshold_selector = ipw.FloatLogSlider(
            value=0.01,
            base=10,
            min=-5,  # max exponent of base
            max=-1,  # min exponent of base
            step=0.1,  # exponent step
            description="Min Energy Diff",
        )

        # Toggle button for filter mode
        self.filter_mode_buttons = ipw.ToggleButtons(
            options={
                label: value
                for label, value in zip(
                    self.FILTER_MODES_DESC, self.FILTER_MODES
                )
            },
            value=self.FILTER_MODES[0],
            description="Filter Mode:",
        )

        # Selected Species (TODO: Make setter/getter)
        self._atomic_number = 2
        self._ion_number = 0

        # Define colors for the transitions based on wavelengths
        self.colorscale = colorscale
        self.cmap = plt.get_cmap(self.colorscale)

        self.compute_level_data()
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "compute_transitions()" explicitly

        # TODO: Revisit later
        self.level_width_scale, self.level_width_offset = 3, 1
        self.transition_width_scale, self.transition_width_offset = 2, 1

    def reset_selected_plot_wavelength_range(self):
        self.min_wavelength = (
            self.line_interaction_analysis[self.filter_mode]
            .lines.loc[self.atomic_number, self.ion_number]
            .wavelength.min()
        )
        self.max_wavelength = (
            self.line_interaction_analysis[self.filter_mode]
            .lines.loc[self.atomic_number, self.ion_number]
            .wavelength.max()
        )
        self.compute_transitions()

    def set_min_plot_wavelength(self, value):
        self.min_wavelength = value
        self.compute_transitions()

    def set_max_plot_wavelength(self, value):
        self.max_wavelength = value
        self.compute_transitions()

    @property
    def max_levels(self):
        return self.max_level_selector.value

    @property
    def level_diff_threshold(self):
        return self.level_diff_threshold_selector.value

    @property
    def filter_mode(self):
        return self.filter_mode_buttons.value

    @property
    def atomic_number(self):
        return self._atomic_number

    @atomic_number.setter
    def atomic_number(self, value):
        self._atomic_number = value
        # TODO: Handle this better
        self.compute_level_data()
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "compute_transitions()" explicitly

    @property
    def ion_number(self):
        return self._ion_number

    @ion_number.setter
    def ion_number(self, value):
        self._ion_number = value
        # TODO: Handle this better
        self.compute_level_data()
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "compute_transitions()" explicitly

    @property
    def atomic_name(self):
        return self.atom_data.loc[self.atomic_number]["name"]

    @property
    def atomic_symbol(self):
        return self.atom_data.loc[self.atomic_number]["symbol"]

    def compute_transitions(self):
        # Get relevant lines for current simulation
        self.line_interaction_analysis[
            self.filter_mode
        ].atomic_number = self.atomic_number
        self.line_interaction_analysis[
            self.filter_mode
        ].ion_number = self.ion_number

        # Get the excitation/de-excitation transitions from LastLineInteraction object
        excit_lines = (
            self.line_interaction_analysis[self.filter_mode]
            .last_line_in.reset_index()
            .groupby(["level_number_lower", "level_number_upper"])
            .agg(
                num_electrons=("line_id", "count"),  # Take count of lines
                wavelength=("wavelength", "first"),  # Take first of wavelengths
            )
            .reset_index()
        )

        deexcit_lines = (
            self.line_interaction_analysis[self.filter_mode]
            .last_line_out.reset_index()
            .groupby(["level_number_lower", "level_number_upper"])
            .agg(
                num_electrons=("line_id", "count"),  # Take count of lines
                wavelength=("wavelength", "first"),  # Take first of wavelengths
            )
            .reset_index()
        )

        # Filter transitions to only include transitions upto the self.max_levels
        excit_lines = excit_lines.loc[
            excit_lines.level_number_upper <= self.max_levels
        ]
        deexcit_lines = deexcit_lines.loc[
            deexcit_lines.level_number_upper <= self.max_levels
        ]

        # Map the levels to merged levels
        excit_lines[
            "merged_level_number_lower"
        ] = excit_lines.level_number_lower.map(self.level_mapping)
        excit_lines[
            "merged_level_number_upper"
        ] = excit_lines.level_number_upper.map(self.level_mapping)
        deexcit_lines[
            "merged_level_number_lower"
        ] = deexcit_lines.level_number_lower.map(self.level_mapping)
        deexcit_lines[
            "merged_level_number_upper"
        ] = deexcit_lines.level_number_upper.map(self.level_mapping)

        # Group by level pairs
        excit_lines = (
            excit_lines.groupby(
                ["merged_level_number_lower", "merged_level_number_upper"]
            )
            .agg(
                wavelength=("wavelength", "mean"),  # Take mean of wavelength
                num_electrons=("num_electrons", "sum"),  # Take sum of counts
            )
            .reset_index()
        )
        deexcit_lines = (
            deexcit_lines.groupby(
                ["merged_level_number_lower", "merged_level_number_upper"]
            )
            .agg(
                wavelength=("wavelength", "mean"),  # Take mean of wavelength
                num_electrons=("num_electrons", "sum"),  # Take sum of counts
            )
            .reset_index()
        )

        # Remove the rows where start and end (merged) level is the same
        excit_lines = excit_lines.loc[
            excit_lines.merged_level_number_lower
            != excit_lines.merged_level_number_upper
        ]
        deexcit_lines = deexcit_lines.loc[
            deexcit_lines.merged_level_number_lower
            != deexcit_lines.merged_level_number_upper
        ]

        # Remove the rows outside the wavelength range for the plot
        excit_lines = excit_lines.loc[
            (excit_lines.wavelength >= self.min_wavelength)
            & (excit_lines.wavelength <= self.max_wavelength)
        ]
        deexcit_lines = deexcit_lines.loc[
            (deexcit_lines.wavelength >= self.min_wavelength)
            & (deexcit_lines.wavelength <= self.max_wavelength)
        ]

        # Compute the standardized log number of electrons for transition line thickness and offset
        excit_log_num_electrons_range = np.log10(
            np.max(excit_lines.num_electrons)
            / np.min(excit_lines.num_electrons)
        )
        excit_lines["standard_log_num_electrons"] = 0
        if excit_log_num_electrons_range > 0:
            excit_lines.standard_log_num_electrons = (
                np.log10(
                    excit_lines.num_electrons
                    / np.min(excit_lines.num_electrons)
                )
                / excit_log_num_electrons_range
            )

        deexcit_log_num_electrons_range = np.log10(
            np.max(deexcit_lines.num_electrons)
            / np.min(deexcit_lines.num_electrons)
        )
        deexcit_lines["standard_log_num_electrons"] = 0
        if deexcit_log_num_electrons_range > 0:
            deexcit_lines.standard_log_num_electrons = (
                np.log10(
                    deexcit_lines.num_electrons
                    / np.min(deexcit_lines.num_electrons)
                )
            ) / deexcit_log_num_electrons_range

        self.excit_lines = excit_lines
        self.deexcit_lines = deexcit_lines

    def compute_level_data(self):
        ### Get energy levels and convert to eV
        raw_energy_levels = self.level_energy_data.loc[
            self.atomic_number, self.ion_number
        ].loc[1 : self.max_levels]

        ### Get level populations
        raw_level_populations = self.level_population_data.loc[
            self.atomic_number, self.ion_number
        ].loc[1 : self.max_levels]

        # Average out the level populations across all zones (TODO: might include an option to select the zone number later)
        raw_level_populations = raw_level_populations.mean(axis=1)

        raw_level_populations = pd.Series(
            raw_level_populations, name="population"
        )

        ### Join level populations and energy values
        raw_level_data = pd.merge(
            raw_energy_levels,
            raw_level_populations,
            left_index=True,
            right_index=True,
        )

        ### Merge the levels if energy difference is too less
        # Get new level numbers
        # TODO: Find a better way to find close levels
        raw_level_data["merged_level_number"] = (
            (raw_level_data["energy"] + 1).pct_change().abs()
            > self.level_diff_threshold
        ).cumsum() + 1

        # Group data with new level numbers
        self.level_data = (
            raw_level_data.reset_index()
            .groupby("merged_level_number")
            .agg(
                energy=(
                    "energy",
                    "mean",
                ),  # Set energy as mean of merged levels
                population=("population", "sum"),
            )
        )  # Add the populations of merged levels

        self.level_mapping = raw_level_data.merged_level_number

        self.merged_max_energy_level = self.level_data.energy.max()

    def display(self):
        ### Create figure and set metadata
        fig = go.FigureWidget()

        fig.update_layout(
            title=f"Grotrian Diagram for {self.atomic_name} {int_to_roman(self.ion_number + 1)}",
            title_x=0.5,
            yaxis_title="Level Number",
            plot_bgcolor="white",
            autosize=False,
            width=700,
            height=700,
            xaxis=dict(showticklabels=False, showgrid=False),
            yaxis=dict(showgrid=False, range=[0, None]),
            showlegend=False,
        )

        ### Create energy level platforms in the figure
        # Create energy tick for ground state separately
        fig.add_annotation(
            x=1.1,
            y=0,
            text=f"{0:.1f} eV",
            showarrow=False,
        )

        # Standardize the level populations to display the widths correctly
        standard_log_population_range = np.log10(
            np.max(self.level_data.population)
            / np.min(self.level_data.population)
        )
        self.level_data["standard_log_population"] = 0
        if standard_log_population_range > 0:
            self.level_data.standard_log_population = (
                np.log10(
                    self.level_data.population
                    / np.min(self.level_data.population)
                )
                / standard_log_population_range
            )
        # Create the energy levels from level data
        for level_number, level_info in self.level_data.iterrows():
            level_width = (
                level_info.standard_log_population * self.level_width_scale
                + self.level_width_offset
            )
            fig.add_trace(
                go.Scatter(
                    x=[0, 1],
                    y=[level_number, level_number],
                    mode="lines+text",
                    hovertemplate=f"Energy: {level_info.energy:.1f} eV<br>"
                    + f"Population: {level_info.population:.2e}"
                    + "<extra></extra>",
                    line=dict(color="black", width=level_width),
                    showlegend=False,
                )
            )

            fig.add_annotation(
                x=1.1,
                y=level_number,
                text=f"{level_info.energy:.1f} eV",
                showarrow=False,
            )

        ### Create transition lines
        wavelength_range = np.log10(self.max_wavelength / self.min_wavelength)
        # Plot excitation transitions
        for _, line_info in self.excit_lines.iterrows():
            lower, upper = (
                line_info.merged_level_number_lower,
                line_info.merged_level_number_upper,
            )
            wavelength, standard_log_num_electrons = (
                line_info.wavelength,
                line_info.standard_log_num_electrons,
            )
            energy_lower, energy_upper = (
                self.level_data.loc[lower].energy,
                self.level_data.loc[upper].energy,
            )

            x_end = (energy_upper - energy_lower) / (
                self.merged_max_energy_level - energy_lower
            )

            color = matplotlib.colors.rgb2hex(
                self.cmap(
                    (np.log10(wavelength / self.min_wavelength))
                    / wavelength_range
                )[:3]
            )

            # Add arrowhead
            fig.add_annotation(
                x=x_end,
                y=upper,  # Start of arrow
                ax=0,
                ay=lower,  # End of arrow
                xref="x",
                yref="y",
                axref="x",
                ayref="y",
                showarrow=True,
                arrowcolor=color,
                arrowhead=2,  # Arrow style
                arrowsize=1,
                arrowwidth=standard_log_num_electrons
                * self.transition_width_scale
                + self.transition_width_offset,  # Adjust width accordingly
            )

        # Plot deexcitation transitions
        for _, line_info in self.deexcit_lines.iterrows():
            lower, upper = (
                line_info.merged_level_number_lower,
                line_info.merged_level_number_upper,
            )
            wavelength, standard_log_num_electrons = (
                line_info.wavelength,
                line_info.standard_log_num_electrons,
            )
            energy_lower, energy_upper = (
                self.level_data.loc[lower].energy,
                self.level_data.loc[upper].energy,
            )

            x_end = (energy_upper - energy_lower) / (
                self.merged_max_energy_level - energy_lower
            )

            color = matplotlib.colors.rgb2hex(
                self.cmap(
                    (np.log10(wavelength / self.min_wavelength))
                    / wavelength_range
                )[:3]
            )

            # Add arrowhead
            fig.add_annotation(
                x=x_end,
                y=lower,  # Start of arrow
                ax=0,
                ay=upper,  # End of arrow
                xref="x",
                yref="y",
                axref="x",
                ayref="y",
                showarrow=True,
                arrowcolor=color,
                arrowhead=2,  # Arrow style
                arrowsize=1,  # Make the arrowhead 2 times bigger
                arrowwidth=standard_log_num_electrons
                * self.transition_width_scale
                + self.transition_width_offset,  # Adjust width accordingly
            )

        # Add a dummy Scatter trace to display colorbar
        tickvalues = np.geomspace(self.min_wavelength, self.max_wavelength, 5)
        ticktext = [f"{val:.1e}" for val in tickvalues]
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    colorscale=self.colorscale,
                    showscale=True,
                    cmin=np.log10(self.min_wavelength),
                    cmax=np.log10(self.max_wavelength),
                    colorbar=dict(
                        title=f"Wavelength ({u.Angstrom})",
                        thickness=5,
                        tickvals=np.log10(tickvalues),
                        ticktext=ticktext,
                        outlinewidth=0,
                        x=1.2,
                    ),
                ),
                hoverinfo="none",
            )
        )

        return fig
