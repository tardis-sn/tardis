import astropy.units as u
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from tardis.util.base import element_symbol2atomic_number
class GrotrianDiagram:
    """
    Prospective plotting interface for GSOC's 2023 idea: Grotrian Diagram plotting. This is the first objective,
    uses Plotly to plot neutral and singly-ionized silicon against velocity. 
    """
    def __init__(self, sim):
        """
        Initializes the GrotrianDiagram object using the simulation object generated from the run_tardis function.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            simulation object generated using the run_tardis function.
        """
        self.model  = sim.model
        self.plasma = sim.plasma
        self.theme_colors = dict(
            light=dict(
                linecolor="#555",
                gridcolor="#e6e6e6",
                zerolinecolor="#fafafa",
                color="#000",
                plot_bgcolor="#fafafa",
                paper_bgcolor="#fafafa",
                title_font_color="#444",
                legendgrouptitle_color="#444",
                bordercolor="#BEC8D9",
                font_color="#000",
            ),
            dark=dict(
                linecolor="#050505",
                gridcolor="#262626",
                zerolinecolor="#000",
                color="#fafafa",
                plot_bgcolor="#000",
                paper_bgcolor="#000",
                title_font_color="#ccc",
                legendgrouptitle_color="#ccc",
                bordercolor="black",
                font_color="#fafafa",
            ),
        )
        #There's only 8 colors here. I don't know how high ionization numbers can get (?)
        self.line_colors = ["#003865", "#005499", "#0071cc","#008dff",
                            "#33a4ff", "#66baff", "#99d1ff", "#007A9B"]
        
        
    def generate_plot(self, atoms_and_ions_list=[("Si",[0, 1])],
                      theme="light", grid_on=False):
        """
        Creates the plot of silicon fractions against velocity.

        Parameters
        ----------
        atoms_and_ions_list : List of pairs (optional)
            first : str
                The name of the atom whose ions are to be plotted
            second : list of int
                The list of different ions from that atom to plot. 
        theme : str, optional
            theme for the plot, by default "light". Could be "dark".
        grid_on : bool, optional
            Whether to show grid or not.
       

        Returns
        -------
        plotly.graph_objs._figure.Figure
            Figure containing the required plots.
            
        Note:
        -----
        For each item in the dictionary above, there will be a different subplot
        
        """
        #TODO: Handle errors when ionic number is too high for element. Or maybe just let pandas do it?
        #TODO: Make legend for subplots align with each using legend_tracegroupgap. --update: won't work, alignment is difficult
        #      Make function return list of figures instead, each will have its own legend.
        #TODO: Add an option to have multiple elements in the same sublot, configure the list and the code accordingly
        

        self.fig = make_subplots(rows=len(atoms_and_ions_list), cols=1)
        self.fig.update_layout(hovermode='x')
            
        
        #Get x-axis values
        v_shells = self.model.velocity.to_value(u.km / u.s)

        if (grid_on):
            zerolinecolor = self.theme_colors[theme]["gridcolor"]
        else:
            zerolinecolor = self.theme_colors[theme]["zerolinecolor"]
            
        self.fig.update_xaxes(
            title="Velocity (km/s)",
            exponentformat="none",
            color=self.theme_colors[theme]["color"],
            linecolor=self.theme_colors[theme]["linecolor"],
            gridcolor=self.theme_colors[theme]["gridcolor"],
            zerolinecolor=zerolinecolor,
            tickangle = 0,
            showgrid=grid_on,
        )
        
        self.fig.update_yaxes(
            title="Ion fraction (unitless)",
            exponentformat="none",
            color=self.theme_colors[theme]["color"],
            linecolor=self.theme_colors[theme]["linecolor"],
            gridcolor=self.theme_colors[theme]["gridcolor"],
            zerolinecolor=zerolinecolor,
            tickformat='.3e',
            tickangle=0,
            showgrid=grid_on,
        )
        
        self.fig.layout.plot_bgcolor = self.theme_colors[theme]["plot_bgcolor"]
        self.fig.layout.paper_bgcolor = self.theme_colors[theme]["paper_bgcolor"]
        
        self.fig.update_layout(
            width=984,
            height=408 * len(atoms_and_ions_list),
            title=dict(
                text='<b>Ionized elements fractions vs. Velocity</b>',
                x=0.5,
                y=0.95,
                font=dict(
                    size=26,
                    color=self.theme_colors[theme]["title_font_color"],
                )
            ),
            legend=dict(
                title_text='<b>(Atom, Charge Number)</b>',
                font=dict(
                    size=16,  
                ),
            ),
            font_color=self.theme_colors[theme]["font_color"],         
        )
        
        #We'll use the "index" variable to keep track of which subplot we're on
        index = 0
        for atom, ions in atoms_and_ions_list:
            atomic_number = element_symbol2atomic_number(atom)
            #Get y-axis values
            fractions = self.get_fractions(self.plasma, atomic_number, ions, len(v_shells))
            for fraction_no in range(len(fractions)):
                self.fig.add_trace(
                        go.Scatter(
                            x=v_shells[:-1],
                            y=fractions[fraction_no],
                            name="<b>(" + atom +", Charge: " +
                                                str(ions[fraction_no]) + ")</b>",
                            mode="markers+lines",
                            showlegend=True,
                            hovertemplate="<b>Velocity</b>: %{x}" +
                                              "<br><b>Fraction</b>: %{y}</br>",
                            line=dict(
                                color=self.line_colors[fraction_no],
                            ),
                        ),
                        col=1,
                        row=index + 1, 
                )
                
            index += 1
   
        return self.fig
    
    
    def get_fractions(self, plasma, atomic_number, ions, shells_no):
        
        """
        Calculates fractions for the y-axis.

        Parameters
        ----------
        plasma: tardis.simulation.Simulation.Plasma
        atomic_number : int
            The atomic number of the atom whose fractions are to be plotted.
        ions : list of int
            The list of different ions from that atom to plot.
        shells_no : int
            Number of velocities to plot against.
        
        Returns
        -------
        list of lists (float)
            Each list is the fractions for that specific ion at different shells/velocities.
        """
        
        #Get sum of densities of an element across same shell but different ions
        #Append fractions of ions/total for each shell
        
        fractions_temp = []
        for shell in range(shells_no - 1):
            total = plasma.ion_number_density[shell][atomic_number][:].sum()
            fractions_temp.append([])
            for ion in ions:
                fractions_temp[shell].append(plasma.ion_number_density[shell][atomic_number][ion]/total)
        
        #Transpose list to go from
        #| ion0_shell0 | ion1_shell0 |             | ion0_shell0 | ion0_shell1 |
        #| ion0_shell1 | ion1_shell1 | etc. to --> | ion1_shell0 | ion1_shell1 | etc.
        
        fractions = np.asarray(fractions_temp).T.tolist()

        return fractions
