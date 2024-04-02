from tardis import run_tardis
from tardis.io.atom_data.util import download_atom_data
from tardis.visualization import SDECPlotter
from tardis.io.model_reader import read_uniform_abundances
from tardis.io.configuration.config_reader import Configuration
import math
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)
import matplotlib.cm as cm
from collections import defaultdict
import plotly.express as px

def abundancesVsVelocity(sim, cmapname = "jet"):
    v_shells = sim.simulation_state.radius.value * 1e-5 / sim.simulation_state.time_explosion.value
    config = Configuration.from_yaml("tardis_example.yml")
    abundances_sections = config.model.abundances
    
    no_shells = len(v_shells)
    abundance_df, isotope_abundance_df = read_uniform_abundances(abundances_sections, no_shells)

    elemental_id = abundance_df[~abundance_df.isna().all(axis=1)].index.tolist()

    isotopic_id = isotope_abundance_df[~isotope_abundance_df.isna().all(axis=1)].index.tolist()

    plot_data = {}
    if elemental_id:
        for i in elemental_id:
            plot_data[atomic_number2element_symbol(i)] = abundance_df.loc[i].tolist()
        if isotopic_id:
            for ls in isotopic_id:
                plot_data[str(atomic_number2element_symbol(ls[0])+" "+int_to_roman(ls[1]))] = isotope_abundance_d.loc[ls].tolist()
    else:
        if isotopic_id:
            for ls in isotopic_id:
                plot_data[str(atomic_number2element_symbol(ls[0])+" "+int_to_roman(ls[1]))] = isotope_abundance_d.loc[ls].tolist()
        else:
            print('No species to print')
            return -1
    #print(plot_data)
    cmap = cm.get_cmap(cmapname, len(plot_data))
    x = v_shells  
    adundance_offset = 0.001
    for idx, (key, values) in enumerate(plot_data.items()):
        color = cmap(idx)
        print(key)
        plt.plot(x, [val+idx*0.001 for val in values], label=key, color=color)

    plt.xlabel('Velocity Shells')
    plt.ylabel('Abundance')
    plt.legend()
    plt.show()
    
def SDECPlot(sim, mode = None):

    plotter = SDECPlotter.from_simulation(sim)
    if mode == "real":
        plotter.generate_plot_mpl(mode).figure.savefig('first_objective_1(real).png')
    elif mode is None:
        plotter.generate_plot_mpl()
    else:
        print("Mode can either be real or empty")
        return -1
    #plotter

def numPacketsLastIntrElements(sim, cmapname = "jet"):
    lines_df = sim.plasma.atomic_data.lines.reset_index().set_index(
            "line_id"
        )
    last_interaction_type = sim.transport.transport_state.last_interaction_type
    last_line_interaction_in_id = sim.transport.transport_state.last_line_interaction_in_id
    last_line_interaction_out_id = sim.transport.transport_state.last_line_interaction_out_id
    

    packets_df = pd.DataFrame(
        {
             "last_interaction_type": last_interaction_type,
             "last_line_interaction_out_id": last_line_interaction_out_id,
             "last_line_interaction_in_id": last_line_interaction_in_id,
        }
    )

    line_mask = (packets_df["last_interaction_type"] > -1) & (
                packets_df["last_line_interaction_in_id"] > -1
            )

    packets_df_line_interaction = packets_df.loc[line_mask].copy()

    packets_df_line_interaction["last_line_interaction_atom"] = (
                lines_df["atomic_number"]
                .iloc[
                    packets_df_line_interaction["last_line_interaction_out_id"]
                ]
                .to_numpy()
            )
        
    packets_df_line_interaction["last_line_interaction_species"] = (
                lines_df["atomic_number"]
                .iloc[
                    packets_df_line_interaction["last_line_interaction_out_id"]
                ]
                .to_numpy()
                * 100
                + lines_df["ion_number"]
                .iloc[
                    packets_df_line_interaction["last_line_interaction_out_id"]
                ]
                .to_numpy()
            )
    plot_data = defaultdict(int)
    no_of_last_interactions = len(packets_df_line_interaction)
    print(packets_df_line_interaction["last_line_interaction_species"].tolist().count(1403))
    for i in range(no_of_last_interactions):
        id_ = packets_df_line_interaction.iloc[i]["last_line_interaction_species"]
        plot_data[id_] += 1
    print(plot_data)
    legend = dict()
    for id_ in plot_data.keys():
        species = ""
        if id_%100 == 0:
            species = atomic_number2element_symbol(int(id_/100))
        else:
            species = atomic_number2element_symbol(int(id_/100))
            species += " " 
            species += int_to_roman((id_%100)+1)
        legend[id_] = species
    final_dict = dict()
    for id_ in plot_data.keys():
        final_dict[legend[id_]] = plot_data[id_]
    plot_df = pd.DataFrame.from_dict(final_dict, orient='index')
    print(plot_df)
    fig = px.bar(plot_df, x=final_dict.keys(), y=final_dict.values(), labels = {'x': 'Species', 'y': 'Number of last interactions'}, color=final_dict.keys())
    fig.show()       
            

download_atom_data('kurucz_cd23_chianti_H_He')

sim = run_tardis("tardis_example.yml", virtual_packet_logging=True)

#SDECPlot(sim)
SDECPlot(sim, "real")
#abundancesVsVelocity(sim)
#numPacketsLastIntrElements(sim)
        


