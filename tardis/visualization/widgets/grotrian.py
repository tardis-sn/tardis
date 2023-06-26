from tardis.analysis import LastLineInteraction
from tardis.util.base import int_to_roman
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
import ipywidgets as ipw

ANGSTROM_SYMBOL = "\u212B"


class GrotrianWidget:
    FILTER_MODES = ("packet_out_nu", "packet_in_nu")
    FILTER_MODES_DESC = ("Emitted Wavelength", "Absorbed Wavelength")

    @classmethod
    def from_simulation(cls, sim, **kwargs):
        """
        Creates a GrotrianWidget object from a Simulation object

        Args:
            sim (tardis.simulation.Simulation): TARDIS simulation object
        """
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
        colorscale="rainbow",
    ):
        """
        Args:
            atom_data (pandas.DataFrame): Mapping from atomic number to symbol and name
            level_energy_data (pandas.Series): Level energies (in eV) indexed by (atomic_number, ion_number, level_number)
            level_population_data (pandas.DataFrame): Level populations indexed by (atomic_number, ion_number, level_number)
                and each column representing the supernova shell
            line_interaction_analysis (tardis.analysis.LastLineInteraction): LastLineInteraction object with the appropriate filters
            colorscale (str, optional): Colorscale for the wavelength info. Defaults to "rainbow".
        """
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
        self._shell = None

        # Define colors for the transitions based on wavelengths
        self.colorscale = colorscale
        self.cmap = plt.get_cmap(self.colorscale)

        self._compute_level_data()
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "_compute_transitions()" explicitly

        # TODO: Revisit later
        self.level_width_scale, self.level_width_offset = 3, 1
        self.transition_width_scale, self.transition_width_offset = 2, 1
        self.arrowhead_size = 9

        self.x_min, self.x_max = 0, 1

    def reset_selected_plot_wavelength_range(self):
        self.min_wavelength = None
        self.max_wavelength = None
        self._compute_transitions()

    def set_min_plot_wavelength(self, value):
        self.min_wavelength = value
        self._compute_transitions()

    def set_max_plot_wavelength(self, value):
        self.max_wavelength = value
        self._compute_transitions()

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

    def set_ion(self, atomic_number, ion_number):
        self._atomic_number = atomic_number
        self._ion_number = ion_number
        # TODO: Handle this better
        self._compute_level_data()
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "_compute_transitions()" explicitly

    @property
    def ion_number(self):
        return self._ion_number

    @property
    def atomic_name(self):
        return self.atom_data.loc[self.atomic_number]["name"]

    @property
    def atomic_symbol(self):
        return self.atom_data.loc[self.atomic_number]["symbol"]

    @property
    def shell(self):
        return self._shell

    def _compute_transitions(self):
        """
        Computes the excitation/de-excitation line transition data for the arrows in the widget
        """
        # Get relevant lines for current simulation
        self.line_interaction_analysis[self.filter_mode].set_ion(
            self.atomic_number, self.ion_number
        )

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

        # Compute default wavelengths if not set by user
        if self.min_wavelength is None:  # Compute default wavelength
            self.min_wavelength = np.min(
                np.concatenate(
                    (excit_lines.wavelength, deexcit_lines.wavelength)
                )
            )
        if self.max_wavelength is None:  # Compute default wavelength
            self.max_wavelength = np.max(
                np.concatenate(
                    (excit_lines.wavelength, deexcit_lines.wavelength)
                )
            )

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
        max_log_num_electrons = np.log10(
            np.max(
                np.concatenate(
                    (excit_lines.num_electrons, deexcit_lines.num_electrons)
                )
            )
        )
        min_log_num_electrons = np.log10(
            np.min(
                np.concatenate(
                    (excit_lines.num_electrons, deexcit_lines.num_electrons)
                )
            )
        )
        log_num_electrons_range = max_log_num_electrons - min_log_num_electrons

        excit_lines["std_log_num_electrons"] = 0
        deexcit_lines["std_log_num_electrons"] = 0

        if log_num_electrons_range > 0:
            excit_lines.std_log_num_electrons = (
                np.log10(excit_lines.num_electrons) - min_log_num_electrons
            ) / log_num_electrons_range
            deexcit_lines.std_log_num_electrons = (
                np.log10(deexcit_lines.num_electrons) - min_log_num_electrons
            ) / log_num_electrons_range

        self.excit_lines = excit_lines
        self.deexcit_lines = deexcit_lines

    def _compute_level_data(self):
        """
        Cleans the level population data for the horizontal platforms in the widget
        """
        ### Get energy levels and convert to eV
        raw_energy_levels = self.level_energy_data.loc[
            self.atomic_number, self.ion_number
        ].loc[0 : self.max_levels]

        ### Get level populations
        raw_level_populations = self.level_population_data.loc[
            self.atomic_number, self.ion_number
        ].loc[0 : self.max_levels]

        # Average out the level populations across all zones, if zone not selected
        if self.shell is None:
            raw_level_populations = raw_level_populations.mean(axis=1)
        else:
            raw_level_populations = raw_level_populations[self.shell]

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
        # TODO: Find a better way to find close levels (less than 0.03 diff in y-coord)
        raw_level_data["merged_level_number"] = (
            (raw_level_data["energy"] + 1).pct_change().abs()
            > self.level_diff_threshold
        ).cumsum()

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

        # Standardize the level populations to display the widths correctly
        std_log_population_range = np.log10(
            self.level_data.population.max()
            / self.level_data.population[self.level_data.population > 0].min()
        )
        self.level_data["std_log_population"] = 0
        if std_log_population_range > 0:
            self.level_data.std_log_population = (
                np.log10(
                    self.level_data.population
                    / self.level_data.population[
                        self.level_data.population > 0
                    ].min()
                )
                / std_log_population_range
            )

        # Create a mapping from original levels to merged levels
        self.level_mapping = raw_level_data.merged_level_number

    def _draw_energy_levels(self):
        """
        Draws the horizontal energy levels on the widget
        """
        y_offset = 0.1  # To offset the energies on the y-axis so as to accomodate energy levele 0 at y=0
        min_log_energy = np.log10(
            self.level_data.energy[self.level_data.energy > 0].min()
        )
        max_log_energy = np.log10(
            self.level_data.energy[self.level_data.energy > 0].max()
        )
        log_energy_range = max_log_energy - min_log_energy

        # Store y-coordinate of levels in the plot (will come in handy for transition arrows as well)
        self.level_data["y_coord"] = (
            np.log10(self.level_data.energy) - min_log_energy
        ) / log_energy_range + y_offset
        # Set y-coordinate of energy level 0 as 0
        self.level_data.y_coord.mask(
            self.level_data.energy == 0, 0, inplace=True
        )

        # Create the energy levels from level data
        for level_number, level_info in self.level_data.iterrows():
            # Add the horizontal line
            self.fig.add_trace(
                go.Scatter(
                    x=np.linspace(self.x_min - 0.05, self.x_max + 0.05, 10),
                    y=level_info.y_coord * np.ones(10),
                    mode="lines",
                    hovertemplate=f"Energy: {level_info.energy:.2e} eV<br>"
                    + f"Population: {level_info.population:.2e}"
                    + "<extra></extra>",
                    line=dict(
                        color="black",
                        width=level_info.std_log_population
                        * self.level_width_scale
                        + self.level_width_offset,
                    )
                    if level_info.population > 0
                    else dict(color="grey", dash="dash"),
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

            # Add label for energy
            self.fig.add_annotation(
                x=self.x_max + 0.15,
                y=level_info.y_coord,
                text=f"n={level_number}",
                showarrow=False,
                xref="x1",
            )

        ### Create width scale
        # Space the populations (log) and corresponding widths (linear) equally
        min_population_idx = self.level_data.population[
            self.level_data.population > 0
        ].idxmin()
        max_population_idx = self.level_data.population.idxmax()

        min_population = self.level_data.population[min_population_idx]
        max_population = self.level_data.population[max_population_idx]

        min_width = (
            self.level_data.std_log_population[min_population_idx]
            * self.level_width_scale
            + self.level_width_offset
        )
        max_width = (
            self.level_data.std_log_population[max_population_idx]
            * self.level_width_scale
            + self.level_width_offset
        )

        scale_granularity = 10
        population_ticks = np.geomspace(
            min_population, max_population, scale_granularity
        )
        width_ticks = np.linspace(min_width, max_width, scale_granularity)
        y_positions = np.linspace(0, 1, scale_granularity)

        # Draw the scale lines
        for population, width, y_pos in zip(
            population_ticks, width_ticks, y_positions
        ):
            self.fig.add_shape(
                type="line",
                line_width=width,
                x0=0.1,
                x1=0.3,
                y0=y_pos,
                y1=y_pos,
                xref="x2",
                yref="y2",
            )
            self.fig.add_annotation(
                x=0.5,
                y=y_pos,
                text=f"{population:.1e}",
                showarrow=False,
                xref="x2",
                yref="y2",
            )
        # Add title of the width bar
        self.fig.add_annotation(
            x=0.4,
            y=-0.05,
            text="Populations",
            showarrow=False,
            xref="x2",
            yref="y2",
        )

        # Add separator
        self.fig.add_shape(
            type="line",
            line=dict(color="grey", dash="dash"),
            line_width=0.5,
            x0=0.7,
            x1=0.7,
            y0=0,
            y1=1,
            xref="x2",
            yref="y2",
        )

    def _draw_transitions(self, is_excitation):
        """
        Draws the transition arrows on the widget
        """
        lines = self.excit_lines if is_excitation else self.deexcit_lines
        wavelength_range = np.log10(self.max_wavelength / self.min_wavelength)
        # Plot excitation transitions
        for _, line_info in lines.iterrows():
            lower, upper = (
                line_info.merged_level_number_lower,
                line_info.merged_level_number_upper,
            )
            wavelength, std_log_num_electrons = (
                line_info.wavelength,
                line_info.std_log_num_electrons,
            )
            energy_lower, energy_upper = (
                self.level_data.loc[lower].energy,
                self.level_data.loc[upper].energy,
            )

            # Get the end x-coordinate (proportional to energy difference between levels)
            merged_max_energy_level = self.level_data.energy.max()
            x_end = (
                (energy_upper - energy_lower)
                * (self.x_max - self.x_min)
                / (merged_max_energy_level - energy_lower)
            )

            # Get the appropriate y-coordinate (computed in _draw_energy_levels)
            y_lower = self.level_data.loc[lower].y_coord
            y_upper = self.level_data.loc[upper].y_coord

            # Get the end arrow color (proportional to log wavelength)
            color = matplotlib.colors.rgb2hex(
                self.cmap(
                    (np.log10(wavelength / self.min_wavelength))
                    / wavelength_range
                )[:3]
            )

            # Draw arrow
            self.fig.add_trace(
                go.Scatter(
                    x=[self.x_min, x_end],
                    y=[y_lower, y_upper]
                    if is_excitation
                    else [y_upper, y_lower],
                    hovertemplate=f"Count: {int(line_info.num_electrons)}<br>"
                    + f"Wavelength: {wavelength:.2e} {ANGSTROM_SYMBOL}"
                    + "<extra></extra>",
                    marker=dict(
                        size=self.arrowhead_size,
                        color=color,
                        symbol="arrow-bar-up",
                        angleref="previous",
                    ),
                    line=dict(
                        color=color,
                        width=std_log_num_electrons
                        * self.transition_width_scale
                        + self.transition_width_offset,
                    ),
                ),
                row=1,
                col=1,
            )

    def _draw_transition_width_scale(self):
        ##### TODO: Correct Draw the arrow width scale lines
        max_num_electrons = np.max(
            np.concatenate(
                (
                    self.excit_lines.num_electrons,
                    self.deexcit_lines.num_electrons,
                )
            )
        )
        min_num_electrons = np.min(
            np.concatenate(
                (
                    self.excit_lines.num_electrons,
                    self.deexcit_lines.num_electrons,
                )
            )
        )

        max_log_num_electrons = np.max(
            np.concatenate(
                (
                    self.excit_lines.std_log_num_electrons,
                    self.deexcit_lines.std_log_num_electrons,
                )
            )
        )
        min_log_num_electrons = np.min(
            np.concatenate(
                (
                    self.excit_lines.std_log_num_electrons,
                    self.deexcit_lines.std_log_num_electrons,
                )
            )
        )

        min_width = (
            min_log_num_electrons * self.transition_width_scale
            + self.transition_width_offset
        )
        max_width = (
            max_log_num_electrons * self.transition_width_scale
            + self.transition_width_offset
        )

        scale_granularity = 10
        num_electrons_ticks = np.geomspace(
            min_num_electrons, max_num_electrons, scale_granularity
        )
        width_ticks = np.linspace(min_width, max_width, scale_granularity)
        y_positions = np.linspace(0, 1, scale_granularity)
        for num_electrons, width, y_pos in zip(
            num_electrons_ticks, width_ticks, y_positions
        ):
            self.fig.add_shape(
                type="line",
                line_width=width,
                x0=0.75,
                x1=0.95,
                y0=y_pos,
                y1=y_pos,
                xref="x2",
                yref="y2",
            )
            self.fig.add_annotation(
                x=1.15,
                y=y_pos,
                text=f"{num_electrons:.1e}",
                showarrow=False,
                xref="x2",
                yref="y2",
            )
        # Add title of the width bar
        self.fig.add_annotation(
            x=0.95,
            y=-0.05,
            text="#Packets",
            showarrow=False,
            xref="x2",
            yref="y2",
        )

    def display(self):
        """
        Parent function to draw the widget (calls other draw methods independently)
        """
        ### Create figure and set metadata
        self.fig = make_subplots(
            rows=1,
            cols=2,
            column_width=[0.6, 0.4],
            specs=[[{}, {}]],
            horizontal_spacing=0.05,
        )

        # Update fig layout
        self.fig.update_layout(
            title=f"Grotrian Diagram for {self.atomic_name} {int_to_roman(self.ion_number + 1)}",
            title_x=0.5,
            plot_bgcolor="white",
            autosize=False,
            width=700,
            height=700,
            xaxis=dict(showticklabels=False, showgrid=False, automargin=True),
            yaxis=dict(title="Energy (eV)", showgrid=False, range=[0, None]),
            margin=dict(),
            showlegend=False,
        )
        self.fig.update_yaxes(
            showticklabels=False, fixedrange=True, row=1, col=2
        )
        self.fig.update_xaxes(
            showticklabels=False, fixedrange=True, row=1, col=2
        )

        ### Create energy level platforms in the figure
        self._draw_energy_levels()

        # Update y-ticks to reflect actual energy instead of log energies
        self.fig.update_xaxes(fixedrange=True, row=1, col=1)
        self.fig.update_yaxes(
            tickmode="array",
            tickvals=self.level_data.y_coord,
            ticktext=[f"{energy:.2e}" for energy in self.level_data.energy],
            fixedrange=True,
            row=1,
            col=1,
        )

        ### Create transition lines
        self._draw_transitions(is_excitation=True)
        self._draw_transitions(is_excitation=False)
        self._draw_transition_width_scale()

        # Add a dummy Scatter trace to display colorbar
        tickvals = np.geomspace(self.min_wavelength, self.max_wavelength, 10)
        ticktext = [f"{val:.1e}" for val in tickvals]
        self.fig.add_trace(
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
                        title=dict(
                            text=f"Wavelength ({ANGSTROM_SYMBOL})", font_size=12
                        ),
                        thickness=5,
                        tickvals=np.log10(tickvals),
                        ticktext=ticktext,
                        outlinewidth=0,
                    ),
                ),
                hoverinfo="none",
            ),
            row=1,
            col=2,
        )

        return self.fig
