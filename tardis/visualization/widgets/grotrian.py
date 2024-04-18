"""
Grotrian Diagram Widget for TARDIS simulation models.

This widget displays a Grotrian Diagram of the last line interactions of the simulation packets
"""
from tardis.analysis import LastLineInteraction
from tardis.util.base import species_tuple_to_string, species_string_to_tuple
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


def is_zero_defined(transform):
    """
    Utility function to decide if a certain transform is defined at zero

    Parameters
    ----------
    transform : function

    Returns
    -------
    bool
        True if transform is defined at 0 else False
    """
    if transform in [np.log, np.log10]:
        return True
    return False


def standardize(
    values,
    transform=lambda x: x,
    min_value=None,
    max_value=None,
    zero_undefined_offset=0,
):
    """
    Utility function to standardize displayed values like wavelengths, num_packets, levels populations to the range [0, 1]
    This helps in computing visual elements like widths, colors, etc.

    Parameters
    ----------
    values : pandas.Series
        The data to standardize
    transform : function, optional
        Transformations like np.log, np.exp, etc. to apply on the data. Defaults to identity
    min_value : float, optional
        The lower bound of the range
    max_value : float, optional
        The upper bound of the range
    zero_undefined_offset : int, optional
        This is useful for log transformation because log(0) is -inf.
        Hence, value=0 gives y=0 while the
        output for other values start at `zero_undefined_offset` (y = log(value) + zero_undefined_offset)
        Default value is 0

    Returns
    -------
    pandas.Series
        Values after standardization
    """
    zero_undefined = is_zero_defined(transform)  # Is function defined at 0?

    if zero_undefined and zero_undefined_offset == 0:
        raise ValueError(
            "If zero of the transformation is undefined, then provide an offset greater than 0"
        )

    # Compute lower and upper bounds of values
    if min_value is None:
        if zero_undefined:
            min_value = (
                values[values > 0].min() if len(values[values > 0]) > 0 else 0
            )
        else:
            min_value = values.min() if len(values) > 0 else 0
    if max_value is None:
        if zero_undefined:
            max_value = (
                values[values > 0].max() if len(values[values > 0]) > 0 else 0
            )
        else:
            max_value = values.max() if len(values) > 0 else 0

    # Apply transformation if given
    transformed_min_value = (
        transform(min_value) if (min_value > 0 or not zero_undefined) else 0
    )
    transformed_max_value = (
        transform(max_value) if (max_value > 0 or not zero_undefined) else 0
    )
    transformed_values = transform(values)

    # Compute range
    value_range = transformed_max_value - transformed_min_value

    # Apply standardization
    if value_range > 0:
        transformed_values = (
            transformed_values - transformed_min_value
        ) / value_range
        if zero_undefined:
            transformed_values = transformed_values + zero_undefined_offset
            transformed_values = np.where(values == 0, 0, transformed_values)
    else:
        # If only single value present in table, then place it at 0
        transformed_values = 0 * values

    return transformed_values


class GrotrianPlot:
    """
    Class for the Grotrian Diagram

    Parameters
    ----------
    atom_data : pandas.DataFrame
        Mapping from atomic number to symbol and name
    level_energy_data : pandas.Series
        Level energies (in eV) indexed by (atomic_number, ion_number, level_number)
    level_population_data : pandas.DataFrame
        Level populations indexed by (atomic_number, ion_number, level_number)
        and each column representing the supernova shell
    line_interaction_analysis : tardis.analysis.LastLineInteraction
        LastLineInteraction object with the appropriate filters

    Configurable Attributes
    -----------------------
    atomic_number : int
        Atomic number of the ion for which the diagram is plotted
        Note: User should set the atomic_number and ion_number together using set_ion function.
    ion_number : int
        Ion number of the ion for which the diagram is plotted
        Note: User should set the atomic_number and ion_number together using set_ion function.
    shell : int or None
        The supernova shell to filter on.
        If None, the level populations are averaged across all shells,
        and all last line interaction are considered
        Default value is None
    max_levels : int
        The maximum number of levels to plot.
        Default value is 10
    level_diff_threshold : float
        The percentage threshold under which levels are merged
        Default value is 1% (0.01)
    min_wavelength : float
        The minimum wavelength allowed for the transitions
    max_wavelength : float
        The maximum wavelength allowed for the transitions
    filter_mode : {"packet_out_nu", "packet_in_nu"}
        The type of wavelength to apply wavelength range filter on
        Default value is packet_out_nu
    y_scale : {"Log", "Linear"}
        The scale to plot the energy levels on the y-axis
        Default value is Log
    cmapname : str
        The name of the colormap used to denote wavelengths. Default value is "rainbow"
    level_width_scale : float
        The multiplier to convert standardized level populations to level widths
        Default value is 3
    level_width_offset : float
        The offset for level widths (to add to the scaled standardized level populations)
        Default value is 1
    transition_width_scale : float
        The multiplier to convert standardized packet count to transition widths
        Default value is 2
    transition_width_offset : float
        The offset for transition widths (to add to the scaled standardized packet counts)
        Default value is 1
    """

    FILTER_MODES = ("packet_out_nu", "packet_in_nu")
    FILTER_MODES_DESC = ("Emitted Wavelength", "Absorbed Wavelength")
    Y_SCALE_OPTION = {"Linear": (lambda x: x), "Log": np.log}

    @classmethod
    def from_simulation(cls, sim, **kwargs):
        """
        Creates a GrotrianPlot object from a Simulation object

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS simulation object

        Returns
        -------
        tardis.visualization.widgets.grotrian.GrotrianPlot
            GrotrianPlot object
        """
        atom_data = sim.plasma.atomic_data.atom_data
        level_energy_data = pd.Series(
            sim.plasma.atomic_data.levels.energy * u.erg.to(u.electronvolt),
            name="energy",
        )
        level_population_data = sim.plasma.level_number_density
        line_interaction_analysis = {
            filter_mode: LastLineInteraction.from_simulation(sim, filter_mode)
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
    ):
        # Set data members
        self._atom_data = atom_data
        self._level_energy_data = level_energy_data
        self._level_population_data = level_population_data
        self._line_interaction_analysis = line_interaction_analysis

        # Max number of levels to display
        self._max_levels = 10

        # Energy difference threshold below which levels are merged
        self._level_diff_threshold = 0.01

        # Filter mode for the wavelength range
        self._min_wavelength = None
        self._max_wavelength = None
        self._filter_mode = self.FILTER_MODES[0]

        # Selected Species
        self._atomic_number = None
        self._ion_number = None
        self._shell = None

        ### Define default parameters for visual elements related to energy levels
        self.level_width_scale, self.level_width_offset = 3, 1
        self._level_width_transform = np.log  # Scale of the level widths
        self._population_spacer = np.geomspace  # To space width bar counts
        ### Scale of the y-axis
        self._y_scale = "Log"
        self._y_coord_transform = self.Y_SCALE_OPTION[self._y_scale]

        ### Define default parameters for visual elements related to transitions
        self.transition_width_scale, self.transition_width_offset = 2, 1
        self._transition_width_transform = np.log  # Scale of the arrow widths
        self._transition_count_spacer = (
            np.geomspace
        )  # To space width bar counts
        self.arrowhead_size = 9

        ### Define default parameters for visual elements related to wavelengths
        self.cmapname = "rainbow"
        self._wavelength_color_transform = np.log  # Scale of wavelength color
        self._wavelength_spacer = np.geomspace  # To space colorbar wavelengths

        # Coordinate end points of levels
        self.x_min, self.x_max = 0, 1

    @property
    def min_wavelength(self):
        return self._min_wavelength

    @min_wavelength.setter
    def min_wavelength(self, value):
        self._min_wavelength = value
        self._compute_transitions()

    @property
    def max_wavelength(self):
        return self._max_wavelength

    @max_wavelength.setter
    def max_wavelength(self, value):
        self._max_wavelength = value
        self._compute_transitions()

    def reset_selected_plot_wavelength_range(self):
        """
        Resets the wavelength range of the selected plot
        """
        self.min_wavelength = None
        self.max_wavelength = None

    @property
    def max_levels(self):
        return self._max_levels

    @max_levels.setter
    def max_levels(self, value):
        assert type(value) is int
        self._max_levels = value
        self._compute_level_data()
        self.reset_selected_plot_wavelength_range()  # calls _compute_transitions() as well

    @property
    def level_diff_threshold(self):
        return self._level_diff_threshold

    @level_diff_threshold.setter
    def level_diff_threshold(self, value):
        assert 0 >= value and value < 1
        self._level_diff_threshold = value
        self._compute_level_data()
        self._compute_transitions()

    @property
    def filter_mode(self):
        return self._filter_mode

    @filter_mode.setter
    def filter_mode(self, value):
        assert value in self.FILTER_MODES

        # Set the atomic_number and ion_number in the appropriate analysis object
        self._line_interaction_analysis[value].set_ion(
            self.atomic_number, self.ion_number
        )
        self._line_interaction_analysis[value].shell = self.shell

        self._filter_mode = value

        self._compute_transitions()

    @property
    def atomic_number(self):
        if self._atomic_number is None:
            raise ValueError("Atomic number is not set")
        return self._atomic_number

    def set_ion(self, atomic_number, ion_number):
        """
        Sets the atomic number and ion number
        """
        assert type(atomic_number) is int and type(ion_number) is int
        if (atomic_number, ion_number) not in self._level_energy_data.index or (
            atomic_number,
            ion_number,
        ) not in self._level_population_data.index:
            raise ValueError(
                "The (atomic_number, ion_number) pair doesn't exist in model"
            )
        self._line_interaction_analysis[self.filter_mode].set_ion(
            atomic_number, ion_number
        )

        self._atomic_number = atomic_number
        self._ion_number = ion_number
        self._compute_level_data()

        # Reset any custom wavelengths if user changes ion
        self.reset_selected_plot_wavelength_range()  # Also computes transition lines so we don't need to call it "_compute_transitions()" explicitly

    @property
    def ion_number(self):
        if self._ion_number is None:
            raise ValueError("Ion number is not set")
        return self._ion_number

    @property
    def atomic_name(self):
        return self._atom_data.loc[self.atomic_number]["name"]

    @property
    def atomic_symbol(self):
        return self._atom_data.loc[self.atomic_number]["symbol"]

    @property
    def shell(self):
        return self._shell

    @shell.setter
    def shell(self, value):
        assert value is None or type(value) is int
        self._line_interaction_analysis[self.filter_mode].shell = value
        self._shell = value
        self._compute_level_data()
        self._compute_transitions()

    @property
    def y_scale(self):
        return self._y_scale

    @y_scale.setter
    def y_scale(self, value):
        assert value in self.Y_SCALE_OPTION
        self._y_scale = value
        self._y_coord_transform = self.Y_SCALE_OPTION[self._y_scale]

    def _compute_transitions(self):
        """
        Computes the excitation/de-excitation line transition data for the arrows in the widget
        """
        ### Get the excitation/de-excitation transitions from LastLineInteraction object
        excite_lines = (
            self._line_interaction_analysis[self.filter_mode]
            .last_line_in.reset_index()
            .groupby(["level_number_lower", "level_number_upper"])
            .agg(
                num_electrons=("line_id", "count"),  # Take count of lines
                wavelength=("wavelength", "first"),  # Take first of wavelengths
            )
            .reset_index()
        )

        deexcite_lines = (
            self._line_interaction_analysis[self.filter_mode]
            .last_line_out.reset_index()
            .groupby(["level_number_lower", "level_number_upper"])
            .agg(
                num_electrons=("line_id", "count"),  # Take count of lines
                wavelength=("wavelength", "first"),  # Take first of wavelengths
            )
            .reset_index()
        )

        ### Filter transitions to only include transitions up to the self.max_levels
        excite_lines = excite_lines.loc[
            excite_lines.level_number_upper <= self.max_levels
        ]
        deexcite_lines = deexcite_lines.loc[
            deexcite_lines.level_number_upper <= self.max_levels
        ]

        ### Map the levels to merged levels
        excite_lines[
            "merged_level_number_lower"
        ] = excite_lines.level_number_lower.map(self.level_mapping)
        excite_lines[
            "merged_level_number_upper"
        ] = excite_lines.level_number_upper.map(self.level_mapping)
        deexcite_lines[
            "merged_level_number_lower"
        ] = deexcite_lines.level_number_lower.map(self.level_mapping)
        deexcite_lines[
            "merged_level_number_upper"
        ] = deexcite_lines.level_number_upper.map(self.level_mapping)

        ### Group by level pairs
        excite_lines = (
            excite_lines.groupby(
                ["merged_level_number_lower", "merged_level_number_upper"]
            )
            .agg(
                wavelength=("wavelength", "mean"),  # Take mean of wavelength
                num_electrons=("num_electrons", "sum"),  # Take sum of counts
            )
            .reset_index()
        )
        deexcite_lines = (
            deexcite_lines.groupby(
                ["merged_level_number_lower", "merged_level_number_upper"]
            )
            .agg(
                wavelength=("wavelength", "mean"),  # Take mean of wavelength
                num_electrons=("num_electrons", "sum"),  # Take sum of counts
            )
            .reset_index()
        )

        ### Remove the rows where start and end (merged) level is the same
        excite_lines = excite_lines.loc[
            excite_lines.merged_level_number_lower
            != excite_lines.merged_level_number_upper
        ]
        deexcite_lines = deexcite_lines.loc[
            deexcite_lines.merged_level_number_lower
            != deexcite_lines.merged_level_number_upper
        ]

        ### Compute default wavelengths if not set by user
        if len(excite_lines) + len(deexcite_lines) > 0:
            if self.min_wavelength is None:  # Compute default wavelength
                self._min_wavelength = np.min(
                    np.concatenate(
                        (excite_lines.wavelength, deexcite_lines.wavelength)
                    )
                )
            if self.max_wavelength is None:  # Compute default wavelength
                self._max_wavelength = np.max(
                    np.concatenate(
                        (excite_lines.wavelength, deexcite_lines.wavelength)
                    )
                )

            ### Remove the rows outside the wavelength range for the plot
            excite_lines = excite_lines.loc[
                (excite_lines.wavelength >= self.min_wavelength)
                & (excite_lines.wavelength <= self.max_wavelength)
            ]
            deexcite_lines = deexcite_lines.loc[
                (deexcite_lines.wavelength >= self.min_wavelength)
                & (deexcite_lines.wavelength <= self.max_wavelength)
            ]

            ### Compute the standardized log number of electrons for arrow line width
            transition_width_coefficient = standardize(
                np.concatenate(
                    (excite_lines.num_electrons, deexcite_lines.num_electrons)
                ),
                transform=self._transition_width_transform,
                zero_undefined_offset=1e-3,
            )
            excite_lines[
                "transition_width_coefficient"
            ] = transition_width_coefficient[: len(excite_lines)]
            deexcite_lines[
                "transition_width_coefficient"
            ] = transition_width_coefficient[len(excite_lines) :]

        self.excite_lines = excite_lines
        self.deexcite_lines = deexcite_lines

    def _compute_level_data(self):
        """
        Computes the level population data for the horizontal platforms in the widget
        """
        ### Get energy levels and convert to eV
        raw_energy_levels = self._level_energy_data.loc[
            self.atomic_number, self.ion_number
        ].loc[0 : self.max_levels]

        ### Get level populations
        raw_level_populations = self._level_population_data.loc[
            self.atomic_number, self.ion_number
        ].loc[0 : self.max_levels]

        ### Average out the level populations across all zones, if zone not selected
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

        ### Merge the levels if energy difference is less than threshold
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

        ### Standardize the level populations to get width coefficient of levels
        self.level_data["level_width_coefficient"] = standardize(
            self.level_data.population,
            transform=self._level_width_transform,
            zero_undefined_offset=1e-3,
        )

        ### Create a mapping from original levels to merged levels
        self.level_mapping = raw_level_data.merged_level_number

    def _draw_energy_levels(self):
        """
        Draws the horizontal energy levels on the widget
        """
        # Transform energies and standardize result to get y-coordinate in range [0, 1]
        self.level_data["y_coord"] = standardize(
            self.level_data.energy,
            transform=self._y_coord_transform,
            zero_undefined_offset=0.1,
        )

        ### Create the energy levels from level data
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
                        width=level_info.level_width_coefficient
                        * self.level_width_scale
                        + self.level_width_offset,
                    )
                    if level_info.population > 0
                    else dict(color="grey", dash="dash"),
                    showlegend=False,
                ),
                row=1,
                col=2,
            )

            # Add label for energy
            self.fig.add_annotation(
                x=self.x_max + 0.1,
                y=level_info.y_coord,
                text=f"{level_number}",
                showarrow=False,
                xref="x2",
                yref="y2",
            )

    def _draw_population_width_scale(self):
        """
        Displays the level population width reference bar
        """
        ### Create width scale
        ### Find lower and upper bounds of populations and corresponding widths
        min_population_idx = self.level_data.population[
            self.level_data.population > 0
        ].idxmin()
        max_population_idx = self.level_data.population.idxmax()

        min_population = self.level_data.population[min_population_idx]
        max_population = self.level_data.population[max_population_idx]

        min_width = (
            self.level_data.level_width_coefficient[min_population_idx]
            * self.level_width_scale
            + self.level_width_offset
        )
        max_width = (
            self.level_data.level_width_coefficient[max_population_idx]
            * self.level_width_scale
            + self.level_width_offset
        )

        ### Space the populations (log) and corresponding widths (linear) equally
        scale_granularity = 10  # Number of scale ticks to display
        population_ticks = self._population_spacer(
            min_population, max_population, scale_granularity
        )
        width_ticks = np.linspace(min_width, max_width, scale_granularity)
        y_positions = np.linspace(0, 1, scale_granularity)

        ### Draw the scale lines
        for population, width, y_pos in zip(
            population_ticks, width_ticks, y_positions
        ):
            self.fig.add_shape(
                type="line",
                line_width=width,
                x0=0.1,
                x1=0.2,
                y0=y_pos,
                y1=y_pos,
                xref="x1",
                yref="y1",
            )
            self.fig.add_annotation(
                x=0.35,
                y=y_pos,
                text=f"{population:.1e}",
                showarrow=False,
                xref="x1",
                yref="y1",
            )
        # Add title of the width bar
        self.fig.add_annotation(
            x=0.28,
            y=-0.08,
            text="Populations",
            showarrow=False,
            xref="x1",
            yref="y1",
        )

    def _draw_transitions(self, is_excitation):
        """
        Draws the transition arrows on the widget
        """
        lines = self.excite_lines if is_excitation else self.deexcite_lines
        lines["color_coefficient"] = standardize(
            lines.wavelength,
            transform=self._wavelength_color_transform,
            zero_undefined_offset=1e-5,
            min_value=self.min_wavelength,
            max_value=self.max_wavelength,
        )

        self._cmap = plt.get_cmap(self.cmapname)  # Float to color map

        ### Plot excitation transitions
        for _, line_info in lines.iterrows():
            lower, upper = (
                line_info.merged_level_number_lower,
                line_info.merged_level_number_upper,
            )
            wavelength, transition_width_coefficient = (
                line_info.wavelength,
                line_info.transition_width_coefficient,
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
            color_coef = line_info.color_coefficient
            color = matplotlib.colors.rgb2hex(self._cmap(color_coef)[:3])

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
                        symbol="arrow-up",
                        angleref="previous",
                    ),
                    line=dict(
                        color=color,
                        width=transition_width_coefficient
                        * self.transition_width_scale
                        + self.transition_width_offset,
                    ),
                ),
                row=1,
                col=2,
            )

    def _draw_transition_width_scale(self):
        """
        Displays the transition count width reference bar
        """
        ### Find lower and upper bounds of num_electrons and corresponding widths
        max_num_electrons = np.max(
            np.concatenate(
                (
                    self.excite_lines.num_electrons,
                    self.deexcite_lines.num_electrons,
                )
            )
        )
        min_num_electrons = np.min(
            np.concatenate(
                (
                    self.excite_lines.num_electrons,
                    self.deexcite_lines.num_electrons,
                )
            )
        )

        max_width_coefficient = np.max(
            np.concatenate(
                (
                    self.excite_lines.transition_width_coefficient,
                    self.deexcite_lines.transition_width_coefficient,
                )
            )
        )
        min_width_coefficient = np.min(
            np.concatenate(
                (
                    self.excite_lines.transition_width_coefficient,
                    self.deexcite_lines.transition_width_coefficient,
                )
            )
        )

        min_width = (
            min_width_coefficient * self.transition_width_scale
            + self.transition_width_offset
        )
        max_width = (
            max_width_coefficient * self.transition_width_scale
            + self.transition_width_offset
        )

        ### Space the num_electrons (log) and corresponding widths (linear) equally
        scale_granularity = 10
        num_electrons_ticks = self._transition_count_spacer(
            min_num_electrons, max_num_electrons, scale_granularity
        )
        width_ticks = np.linspace(min_width, max_width, scale_granularity)
        y_positions = np.linspace(0, 1, scale_granularity)

        ### Draw the width bar
        for num_electrons, width, y_pos in zip(
            num_electrons_ticks, width_ticks, y_positions
        ):
            self.fig.add_shape(
                type="line",
                line_width=width,
                x0=0.65,
                x1=0.75,
                y0=y_pos,
                y1=y_pos,
                xref="x1",
                yref="y1",
            )
            self.fig.add_annotation(
                x=0.9,
                y=y_pos,
                text=f"{num_electrons:.1e}",
                showarrow=False,
                xref="x1",
                yref="y1",
            )
        # Add title of the width bar
        self.fig.add_annotation(
            x=0.83,
            y=-0.08,
            text="#Packets",
            showarrow=False,
            xref="x1",
            yref="y1",
        )

    def _draw_transition_color_scale(self):
        """
        Displays the transition wavelength colorbar
        """
        # Add a dummy Scatter trace to display colorbar
        tickvals = self._wavelength_spacer(
            self.min_wavelength, self.max_wavelength, 10
        )
        ticktext = [f"{val:.1e}" for val in tickvals]
        self.fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    colorscale=self.cmapname,
                    showscale=True,
                    cmin=self._wavelength_color_transform(self.min_wavelength),
                    cmax=self._wavelength_color_transform(self.max_wavelength),
                    colorbar=dict(
                        title=dict(
                            text=f"Wavelength ({ANGSTROM_SYMBOL})<br>&nbsp;",
                            font_size=12,
                        ),
                        thickness=5,
                        tickvals=self._wavelength_color_transform(tickvals),
                        ticktext=ticktext,
                        outlinewidth=0,
                    ),
                ),
                hoverinfo="none",
            ),
            row=1,
            col=2,
        )

    def display(self):
        """
        Function to draw the plot and the reference scales (calls other draw methods independently)
        """
        ### Create figure and set metadata
        self.fig = go.FigureWidget(
            make_subplots(
                rows=1,
                cols=2,
                column_width=[0.3, 0.7],
                specs=[[{}, {}]],
                horizontal_spacing=0.14,
            )
        )

        # Update fig layout
        self.fig.update_layout(
            title=(
                f"Energy Level Diagram for {self.atomic_name} {int_to_roman(self.ion_number + 1)} "
                f"(Shell: {self.shell if self.shell is not None else 'All'})"
            ),
            title_x=0.5,
            plot_bgcolor="white",
            autosize=False,
            width=1000,
            height=700,
            margin=dict(),
            showlegend=False,
        )

        # Remove ticklabels in the reference bars subplot
        self.fig.update_yaxes(
            showticklabels=False, fixedrange=True, row=1, col=1
        )
        self.fig.update_xaxes(
            showticklabels=False, fixedrange=True, row=1, col=1
        )

        ### Create energy level platforms and width reference scale
        self._draw_energy_levels()
        self._draw_population_width_scale()

        # Remove ticklabels from x-axis
        self.fig.update_xaxes(
            showticklabels=False, fixedrange=True, row=1, col=2
        )
        # Update y-ticks to reflect actual energy values
        self.fig.update_yaxes(
            title=dict(text="Energy (eV)", standoff=5),
            range=[0, None],
            tickmode="array",
            tickvals=self.level_data.y_coord,
            ticktext=[f"{energy:.2e}" for energy in self.level_data.energy],
            fixedrange=True,
            row=1,
            col=2,
        )

        # Add separator between width scales
        self.fig.add_shape(
            type="line",
            line=dict(color="grey", dash="dash"),
            line_width=0.5,
            x0=0.55,
            x1=0.55,
            y0=0,
            y1=1,
            xref="x1",
            yref="y1",
        )

        ### Create transition lines and corresponding width and color scales
        if len(self.excite_lines) > 0:
            self._draw_transitions(is_excitation=True)

        if len(self.deexcite_lines) > 0:
            self._draw_transitions(is_excitation=False)

        if len(self.excite_lines) + len(self.deexcite_lines) > 0:
            self._draw_transition_width_scale()
            self._draw_transition_color_scale()

        return self.fig


class GrotrianWidget:
    """
    A wrapper class for the Grotrian Diagram, containing the Grotrian Plot and the IpyWidgets

    Parameters
    ----------
    plot : tardis.visualization.widgets.grotrian.GrotrianPlot
        GrotrianPlot object
    num_shells : int
        Number of shells in the sim.simulation_state.v_inner
    """

    @classmethod
    def from_simulation(cls, sim, **kwargs):
        """
        Creates a GrotrianWidget object from a Simulation object

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS simulation object

        Returns
        -------
        tardis.visualization.widgets.grotrian.GrotrianWidget
            GrotrianWidget object
        """
        plot = GrotrianPlot.from_simulation(sim, **kwargs)
        num_shells = len(sim.simulation_state.v_inner)
        return cls(plot, num_shells, **kwargs)

    def __init__(self, plot, num_shells, **kwargs):
        self.plot = plot
        self.num_shells = num_shells

        species_list = self._get_species()
        self.ion_selector = ipw.Dropdown(
            options=species_list,
            index=0,
            description="Ion",
        )
        self.plot.set_ion(*species_string_to_tuple(self.ion_selector.value))
        self.ion_selector.observe(
            self._ion_change_handler,
            names="value",
        )
        self.ion_selector.observe(
            self._wavelength_resetter,
            names="value",
        )

        shell_list = ["All"] + [str(i) for i in range(1, num_shells + 1)]
        self.shell_selector = ipw.Dropdown(
            options=shell_list,
            index=0,
            description="Shell",
        )
        self.shell_selector.observe(
            lambda change: self._change_handler(
                "shell", None if change["new"] == "All" else int(change["new"])
            ),
            names="value",
        )
        self.shell_selector.observe(
            self._wavelength_resetter,
            names="value",
        )

        self.max_level_selector = ipw.BoundedIntText(
            value=plot.max_levels,
            min=1,
            max=40,
            step=1,
            description="Max Levels",
        )
        self.max_level_selector.observe(
            lambda change: self._change_handler("max_levels", change["new"]),
            names="value",
        )
        self.max_level_selector.observe(
            self._wavelength_resetter,
            names="value",
        )

        self.y_scale_selector = ipw.ToggleButtons(
            options=GrotrianPlot.Y_SCALE_OPTION.keys(),
            index=1,
            description="Y-Scale",
            layout=ipw.Layout(width="auto"),
            style={"button_width": "100px"},
        )
        self.y_scale_selector.observe(
            lambda change: self._change_handler("y_scale", change["new"]),
            names="value",
        )

        self.wavelength_range_selector = ipw.FloatRangeSlider(
            value=[self.plot.min_wavelength, self.plot.max_wavelength],
            min=self.plot.min_wavelength,
            max=self.plot.max_wavelength,
            step=0.1,
            description="Wavelength",
            layout=ipw.Layout(width="605px"),
            readout_format=".1e",
        )
        self.wavelength_range_selector.observe(
            self._wavelength_change_handler,
            names="value",
        )

    def _get_species(self):
        """
        Computes the ions list for the ion dropdown of the plot
        """
        line_interaction_analysis = self.plot._line_interaction_analysis
        selected_species_group = line_interaction_analysis[
            self.plot.filter_mode
        ].last_line_in.groupby(["atomic_number", "ion_number"])

        if selected_species_group.groups:
            selected_species_symbols = [
                species_tuple_to_string(item)
                for item in selected_species_group.groups.keys()
            ]
        return selected_species_symbols

    def _change_handler(self, attribute, value):
        """
        Generic function to update the configurable attributes of GrotrianPlot object

        Parameters
        ----------
        attribute : str
            The name of the attribute of the GrotrianPlot object
        value :
            The new value of the attribute
        """
        index = self.fig.children.index(self.plot.fig)
        setattr(self.plot, attribute, value)  # Set the value of the attribute

        # Set the updated plot in the figure
        children_list = list(self.fig.children)
        children_list[index] = self.plot.display()
        self.fig.children = tuple(children_list)

    def _ion_change_handler(self, change):
        """
        Function to update ion of GrotrianPlot object

        Parameters
        ----------
        change : dict
            Change information of the event
        """
        atomic_number, ion_number = species_string_to_tuple(change["new"])
        index = self.fig.children.index(self.plot.fig)
        self.plot.set_ion(atomic_number, ion_number)

        # Set the updated plot in the figure
        children_list = list(self.fig.children)
        children_list[index] = self.plot.display()
        self.fig.children = tuple(children_list)
        # self._wavelength_resetter()

    def _wavelength_change_handler(self, change):
        """
        Function to update the wavelength range of GrotrianPlot object

        Parameters
        ----------
        change : dict
            Change information of the event
        """
        min_wavelength, max_wavelength = change["new"]
        index = self.fig.children.index(self.plot.fig)
        setattr(self.plot, "min_wavelength", min_wavelength)
        setattr(self.plot, "max_wavelength", max_wavelength + 1)

        # Set the updated plot in the figure
        children_list = list(self.fig.children)
        children_list[index] = self.plot.display()
        self.fig.children = tuple(children_list)

    def _wavelength_resetter(self, change):
        """
        Resets the range of the wavelength slider whenever the ion, level or shell changes
        """
        min_wavelength = self.plot.min_wavelength
        max_wavelength = self.plot.max_wavelength

        if min_wavelength is None or max_wavelength is None:
            self.wavelength_range_selector.layout.visibility = "hidden"
            return

        elif min_wavelength == max_wavelength:
            self.wavelength_range_selector.layout.visibility = "visible"
            self.wavelength_range_selector.disabled = True
        else:
            self.wavelength_range_selector.layout.visibility = "visible"
            self.wavelength_range_selector.disabled = False

        self.wavelength_range_selector.min = 0.0
        self.wavelength_range_selector.max = max_wavelength
        self.wavelength_range_selector.min = min_wavelength
        self.wavelength_range_selector.value = [
            self.wavelength_range_selector.min,
            self.wavelength_range_selector.max,
        ]

    def display(self):
        """
        Function to render the Grotrian Widget containing the plot and IpyWidgets together
        """
        fig = self.plot.display()
        self.fig = ipw.VBox(
            [
                ipw.HBox(
                    [
                        self.ion_selector,
                        self.shell_selector,
                        self.max_level_selector,
                    ]
                ),
                ipw.HBox(
                    [self.y_scale_selector, self.wavelength_range_selector]
                ),
                fig,
            ]
        )
        return self.fig
