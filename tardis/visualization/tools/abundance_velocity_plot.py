"""
Abundance vs Velocity Plot for TARDIS simulation models.

This plot is the first objective for Tardis Custom Abundance Widget project idea.
There are two function for plotting this graph using different libraries.
"""

import pandas
import numpy as np
import matplotlib.pyplot as plt
from tardis.util.base import atomic_number2element_symbol
import plotly.graph_objects as go
import numpy as np


class AVData:
    """The data of simulation model used by Abundace Velocity (AV) Plot.

    This preprocesses the data required by AVPlotter class for doing
    calculations and plotting
    """

    def __init__(self, av_df, velocity, elements):

        """
        Initialize the AVData with required properties of simulation model.

        Parameters
        ----------
        av_df : pd.DataFrame
            Data about the atomic abundance present in simulation model
        velocity : numpy.ndarray
             An np array that stores inner shell velocities  in cm /s
        elements : numpy.ndarray
             It stores values of all elements present in the model
        """
        self.av_df = av_df
        self.velocity = velocity
        self.elements = elements

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of AVData from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation


        Returns
        -------
        AVData
        """
        av_df = (sim.model.abundance).T
        velocity = sim.model.velocity[:-1]
        elements = av_df.columns.values
        return cls(av_df=av_df, velocity=velocity, elements=elements)


class AVPlotter:
    """
    Plotting interface for Abundance vs Velocity (AV) Plot.

    It performs necessary calculations to generate AV Plot for a simulation
    model, and allows to plot it in matplotlib and plotly.
    """

    def __init__(self, data):
        """
        Initialize the AVPlotter with required data of simulation model.

        Parameters
        ----------
        data : DataFrame
            Abundance and velocity data required for AV plot
        """
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of AVPlotter from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        AVPlotter
        """
        return cls(data=AVData.from_simulation(sim))

    def abundance_velocity_mpl(self):
        """
        Plots the abundance vs velocity of elements provided simulation result using matplotlib

        Returns
        -------
        None
        """

        plt.figure(figsize=(20, 6))

        # Plot line graph for each element
        for atomic_number in self.data.elements:
            plt.plot(
                self.data.velocity,
                self.data.av_df[atomic_number],
                marker="o",
                markersize=12,
                linewidth=4,
                label=atomic_number2element_symbol(atomic_number),
            )

        # Edit the layout
        plt.title("Abundance vs Velocity", fontsize=18)
        plt.xlabel("Velocity ( Inner Shell Boundary in cm/s)")
        plt.ylabel("Abundance (Fraction)")
        plt.legend(loc="best", title="Atomic number", fontsize=10)
        plt.show()

    def abundance_velocity_plotly(self):
        """
        Plots the abundance vs velocity of elements provided simulation result using plotly.

        Returns
        -------
        plotly.graph_objs._figure.Figure
            Figure object on which AV Plot is created
        """
        fig = go.Figure()

        # Plot line graph for each element
        for atomic_number in self.data.elements:
            fig.add_trace(
                go.Scatter(
                    x=self.data.velocity,
                    y=self.data.av_df[atomic_number],
                    mode="lines+markers",
                    name=atomic_number2element_symbol(atomic_number),
                )
            )

        # Edit the layout
        fig.update_layout(
            title="Abundance vs Velocity",
            xaxis_title="Velocity ( Inner Shell Boundary in cm/s )",
            yaxis_title="Abundance (Fraction)",
        )
        fig.update_xaxes(tickformat="e")
        fig.show()

        return fig
