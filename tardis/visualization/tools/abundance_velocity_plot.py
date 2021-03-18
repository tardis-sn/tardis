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
    

def abundance_velocity_matplot(sim):
    """
    Plots the abundance vs velocity of elements provided simulation result using matplotlib.
    Parameters
    ----------
    sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
    Returns
    -------
    None
    """
    plt.figure(figsize=(20,6))

    #Data for plot
    plotdf = (sim.model.abundance).T
    velocity = sim.model.velocity
    elements=plotdf.columns.values

    #Plot line graph for each element
    for atomic_number in elements:
        plt.plot( velocity[:-1], plotdf[atomic_number],marker='o', markersize=12, linewidth=4,label=atomic_number2element_symbol(atomic_number))
    
    # Edit the layout
    plt.title("Abundance vs Velocity", fontsize=18)
    plt.xlabel("Velocity ( Inner Shell Boundary in " + str(sim.model.velocity.unit)+")")
    plt.ylabel("Abundance (Fraction)")
    plt.legend(loc="best",title="Atomic number", fontsize=10)
    plt.show()


def abundance_velocity_plotly(sim): 
    """
    Plots the abundance vs velocity of elements provided simulation result using plotly.
    Parameters
    ----------
    sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
    Returns
    -------
    None
    """
    fig = go.Figure()

    #Data for plot
    plotdf = (sim.model.abundance).T
    velocity = sim.model.velocity
    elements=plotdf.columns.values

    #Plot line graph for each element
    for atomic_number in elements:
        fig.add_trace(go.Scatter(x=velocity, y=plotdf[atomic_number],
                    mode='lines+markers',
                    name=atomic_number2element_symbol(atomic_number)))

    # Edit the layout
    fig.update_layout(title="Abundance vs Velocity",
                   xaxis_title="Velocity ( Inner Shell Boundary in " + str(sim.model.velocity.unit)+")",
                   yaxis_title="Abundance (Fraction)"
                   )
    fig.update_xaxes(tickformat='e') 
    fig.show()


