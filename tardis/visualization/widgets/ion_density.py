"""Class to create and display Ion Density Widget."""
from astropy import units as u
import numpy as np
from plotly import graph_objects as go
import ipywidgets as ipw

from tardis.util.base import (
    atomic_number2element_symbol,
    species_tuple_to_string,
)


class IonDensityWidget:
    """
    Widget to plot the fractional ion density of each shell as a stepped plot

    It allows selection of the element/isotope through a dropdown,
    and the plot is displayed for the selected element/isotope.
    The user can also show/hide the plotlines of each ion by clicking on legend.
    """

    def __init__(
        self,
        ion_density,
        velocity,
    ):
        """
        Initialize the IonDensityWidget with ion density information and shell velocities.

        Parameters
        ----------
        ion_density : pd.DataFrame
            Data about the ion densities present in simulation model's plasma
        velocity : astropy.Quantity
            Velocities at the boundaries of each shell
            This has the inner_boundary velocity at index zero and
            the outer_boundary velocity at the last index
        """
        self.velocity = velocity.to(u.km / u.s)
        self.ion_density = ion_density

        self.atomic_numbers = list(
            ion_density.index.get_level_values("atomic_number").unique()
        )
        self.atomic_symbols = list(
            map(atomic_number2element_symbol, self.atomic_numbers)
        )

        # Widgets ------------------------------------------------
        # Default value is Silicon
        si_atomic_number = 14
        si_symbol = atomic_number2element_symbol(si_atomic_number)
        si_index = self.atomic_symbols.index(si_symbol)

        self.element_dropdown = ipw.Dropdown(
            options=self.atomic_symbols,
            value=si_symbol,
            description="Selected Element",
            index=si_index,
        )

        self.step_plot = self._create_step_plot()
        self._update_step_plot(si_atomic_number)

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of IonDensityWidget from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        IonDensityWidget object
        """
        return cls(
            ion_density=sim.plasma.ion_number_density,
            velocity=sim.model.velocity,
        )

    def _create_step_plot(self):
        """
        Creates the skeleton of the Ion Density plot (without the plot lines)
        It defines the axis names, title, and layout

        Returns
        -------
        plotly.graph_objects.FigureWidget
        """
        fig = go.FigureWidget()

        # Update fig layout and metadata
        fig.update_yaxes(
            type="log", exponentformat="power", title="Ion Density"
        )
        fig.update_xaxes(
            exponentformat="none", title=f"Velocity ({self.velocity.unit})"
        )
        fig.update_layout(
            title_text="Fractional Ion Density vs Velocity",
            title_x=0.5,
            legend=dict(y=0.5, font_size=16),
        )
        return fig

    def _update_step_plot(self, atomic_number):
        """
        Draws the plot lines for the selected element in the dropdown
        It filters out the ions which have zero densities in all shells
        """
        # Filter densities by atomic number
        ion_density = self.ion_density.loc[atomic_number]

        # Get total ion density for the selected atomic number
        atomic_density = ion_density.sum(axis=0)
        # Duplicate innermost density of first shell to represent photosphere
        atomic_density = atomic_density[:1].append(
            atomic_density, ignore_index=True
        )

        # Get the list of ion numbers
        ion_numbers = list(
            ion_density.index.get_level_values("ion_number").unique()
        )

        plotlines = []
        # Add plots for each ion
        for ion_number in ion_numbers:
            # Get ion name (with Roman numerals)
            ion_name = species_tuple_to_string((atomic_number, ion_number))

            # Get the shell densities of this ion
            shell_density = ion_density.loc[ion_number]
            # Duplicate innermost density of first shell to represent photosphere
            shell_density = shell_density[:1].append(
                shell_density, ignore_index=True
            )

            # Skip adding trace if all shell densities are zero
            if np.all(shell_density == 0):
                continue

            # Add plot line
            plotlines.append(
                go.Scatter(
                    x=self.velocity,
                    y=shell_density / atomic_density,
                    name=ion_name,
                    line_shape="vh",
                    hovertemplate=ion_name
                    + "<br>Velocity: %{x:.0f}<br>Density: %{y}",
                )
            )

        self.step_plot.data = []
        self.step_plot.add_traces(plotlines)
        self.step_plot.update_traces(mode="lines+markers")
        self.step_plot.update_layout(
            title_text=f"{atomic_number2element_symbol(atomic_number)} Fractional Ion Density vs Velocity"
        )

    def _element_dropdown_handler(self, change):
        """
        Event handler for selection in element_dropdown.

        This method has the expected signature of the callback function
        passed to :code:`observe` method of ipywidgets as explained in
        `their docs <https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html#Signatures>`_.
        """
        self._update_step_plot(
            self.atomic_numbers[self.element_dropdown.index],
        )

    def display(self):
        """
        Display the fully-functional ion density widget.

        It puts together all component widgets nicely together and enables
        interaction between all the components.

        Returns
        -------
        ipywidgets.Box
            Ion density widget containing all component widgets
        """
        # Set widths of widgets
        self.element_dropdown.style.description_width = "initial"
        self.element_dropdown.style.font_size = "14pt"
        # self.element_dropdown.layout.width = "90%"

        # Attach event listeners to widgets
        self.element_dropdown.observe(
            self._element_dropdown_handler, names="index"
        )

        box_layout = ipw.Layout(
            display="flex",
            flex_flow="column",
            align_items="center",
            margin="10px 0px 0px 5px",
        )

        return ipw.VBox(
            [
                self.element_dropdown,
                self.step_plot,
            ],
            layout=box_layout,
        )
