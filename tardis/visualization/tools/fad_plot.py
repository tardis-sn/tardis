"""
Fractional AbunDance vs velocity for different shells (FAD) Plot for TARDIS simulation models.
"""
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from tardis.util.base import atomic_number2element_symbol


class FADPlotter:
    """
    Plotting interface for Fractional AbunDance vs velocity for different shells (FAD) Plot.

    It offers functions to generate FAD Plot for a simulation
    model, and allows to plot it in matplotlib and plotly.
    """

    def __init__(self, data):
        """
        Initialize the FADPlotter with required data of simulation model.

        Parameters
        ----------
        data : DataFrame
            Abundance and velocity data required for FAD plot
        """
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of FAD from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        ----------
        FADPlotter
        """
        data = sim.model.abundance.T
        data.index = (
            sim.model.velocity[1:] / 100000
        )  # inplace the columns with velocity of shell outer boundary
        data.columns = data.columns.map(
            atomic_number2element_symbol
        )  # map atomic number to element symbol

        return cls(data)

    def generate_plot_mpl(
        self,
        figsize=(15, 5),
    ):
        """
        Generate Spectral element DEComposition (SDEC) Plot using matplotlib.

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
        """

        self.data.plot(figsize=figsize, marker="o", fontsize=12)
        plt.title("Abundance as a function of velocity", size=16)
        plt.ylabel("Fractional Abundance", fontsize=12)
        plt.xlabel("Velocity of Shell Outer Boundary (km/s)", fontsize=12)
        plt.legend(loc=1, title="Element", fontsize=12)

        return plt.gca()

    def generate_plot_ply(self):
        """
        Generate interactive Spectral element DEComposition (SDEC) Plot using plotly.

        Returns
        -------
        plotly.graph_objs._figure.Figure
        """
        title = "Fractional abundance vs velocity of different shells"
        fig = go.Figure()

        for col in self.data.columns:
            fig.add_trace(
                go.Line(
                    x=self.data.index,
                    y=self.data[col],
                    name=self.data[col].name,
                )
            )

        fig.update_layout(
            xaxis=dict(
                title="Velocity of Shell Outer Boundary (km/s)",
                exponentformat="none",
            ),
            yaxis=dict(title="Fractional Abundance", exponentformat="e"),
            height=500,
            title=title,
        )

        return fig
