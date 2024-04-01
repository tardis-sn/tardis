from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.io.atom_data.util import download_atom_data
from collections import defaultdict
from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)
import matplotlib.cm as cm
from matplotlib import colormaps
import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np
import random

def parse_species_list(species_list):
    plot_species_info = []
    if species_list is not None:
        if any(char.isdigit() for char in " ".join(species_list)):
            raise ValueError("All species must be in Roman numeral form, e.g. Si II")
        
        full_species_list = []
        for species in species_list:
            if "-" in species:
                element = species.split(" ")[0]
                first_ion_numeral = roman_to_int(species.split(" ")[-1].split("-")[0])
                second_ion_numeral = roman_to_int(species.split(" ")[-1].split("-")[-1])
                for ion_number in np.arange(first_ion_numeral, second_ion_numeral + 1):
                    full_species_list.append(f"{element} {int_to_roman(ion_number)}")
            else:
                full_species_list.append(species)
        
        requested_species_ids = []
        keep_colour = []
        
        for species in full_species_list:
            if " " in species:
                requested_species_ids.append([species_string_to_tuple(species)[0] * 100 +
                                              species_string_to_tuple(species)[1]])
                plot_species_info.append(requested_species_ids[-1][0])
            else:
                atomic_number = element_symbol2atomic_number(species)
                requested_species_ids.append([atomic_number * 100 + ion_number 
                                              for ion_number in np.arange(atomic_number)])
                keep_colour.append(atomic_number)
                plot_species_info.append(atomic_number*100)
        
        requested_species_ids = [species_id for temp_list in requested_species_ids 
                                 for species_id in temp_list]
        
        return requested_species_ids, keep_colour, plot_species_info
    
    else:
        return None, None

download_atom_data('kurucz_cd23_chianti_H_He')
'''
config = Configuration.from_yaml("tardis_example.yml")


config["montecarlo"]["tracking"]["track_rpacket"]=True
config["montecarlo"]["seed"]= 457
config["montecarlo"]["no_of_packets"]=10
config["montecarlo"]["iterations"]=1
config["montecarlo"]["last_no_of_packets"]=15
config["montecarlo"]["no_of_virtual_packets"]=3

sim = run_tardis(config)
'''
sim = run_tardis("tardis_example.yml", virtual_packet_logging=True)

lines_df = sim.plasma.atomic_data.lines.reset_index().set_index(
            "line_id"
        )
last_interaction_type = sim.transport.transport_state.last_interaction_type
last_line_interaction_in_id = sim.transport.transport_state.last_line_interaction_in_id
last_line_interaction_out_id = sim.transport.transport_state.last_line_interaction_out_id
last_shell = sim.transport.transport_state.last_line_interaction_shell_id
time_explosion = sim.simulation_state.time_explosion.value
v_shells = sim.simulation_state.radius.value * 1e-5 / sim.simulation_state.time_explosion.value

packets_df = pd.DataFrame(
    {
         "last_interaction_type": last_interaction_type,
         "last_line_interaction_out_id": last_line_interaction_out_id,
         "last_line_interaction_in_id": last_line_interaction_in_id,
         "last_shell": last_shell
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

last_interaction_velocity = []
for i in packets_df_line_interaction["last_shell"]:
    last_interaction_velocity.append(v_shells[i])

packets_df_line_interaction["last_interaction_velocity"] = last_interaction_velocity

#print(packets_df_line_interaction)


plot_data = defaultdict(lambda : defaultdict(int))

no_of_interactions = len(packets_df_line_interaction)

for i in range(no_of_interactions):
    plot_data[int(packets_df_line_interaction.iloc[i]["last_line_interaction_species"])][packets_df_line_interaction.iloc[i]["last_interaction_velocity"]] += 1
    
#print(plot_data)

species_list = ["Si III", "Si II", "Si I", "Mg II", "O"]

no_species_required = len(species_list)

species_id = parse_species_list(species_list)

required_plot_data = defaultdict(lambda : defaultdict(int))

for i in species_id[0]:
    if plot_data[i]:
        required_plot_data[i] = plot_data[i]

merged_dict = defaultdict(lambda : defaultdict(int))
ind = []
for key in required_plot_data.keys():
    if int(key/100) in species_id[1]:
        new_key = int(key/100)*100
        for keys in required_plot_data[key].keys():
            merged_dict[new_key][keys] += required_plot_data[key][keys]
        ind.append(key)
for i in ind:
    del required_plot_data[i]
for key in merged_dict.keys():
    required_plot_data[key] = merged_dict[key]
       
print(required_plot_data)      

legend = dict()

for i in range(no_species_required):
    legend[species_list[i]] = species_id[2][i]
 
#print(legend)   
cmapname = "jet"
cmap = cm.get_cmap(cmapname, len(species_list))
    
#print(cmap)
'''
rounded_dict = defaultdict(lambda : defaultdict(int))

for key1, inner_dict in required_plot_data.items():
    for inner_key, value in inner_dict.items():
        rounded_key = round(inner_key, -3)
        rounded_dict[key1][rounded_key] += value

#print(rounded_dict)
'''
x_offset = 50

for i, species in enumerate(species_list):
    x_values = [x + i*x_offset for x in required_plot_data[legend[species]].keys()] 
    y_values = list(required_plot_data[legend[species]].values())
    color = cmap(i) 
    plt.vlines(x_values, ymin=0, ymax=y_values, label=species, color=color)

plt.xlabel('Last Interaction Velocity (km/s)')
plt.ylabel('Packet Count')

plt.legend()

plt.savefig('plot.png')
