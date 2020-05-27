from tardis.base import run_tardis
from tardis.io.atom_data.util import download_atom_data
from tardis.util.base import atomic_number2element_symbol, species_tuple_to_string

import pandas as pd
import numpy as np
import qgrid
import ipywidgets as ipw

# TODO: Fetch saved Simulation object from memory
sim = run_tardis('tardis_example.yml')  # TODO: Make the file path work

class ShellInfo():
    def __init__(self, sim_model=sim):
        self.sim_model = sim_model

    def shells_data(self):
        shells_temp_w = pd.DataFrame({'Rad. Temp.': self.sim_model.model.t_rad,
                                'W': self.sim_model.model.w},
                                index=range(1, 21))
        shells_temp_w.index.name = 'Shell No.'
         # Format to string to make qgrid show values in scientific notations
        return shells_temp_w.applymap(lambda x: '{:.6e}'.format(x))

    def Z_count(self, shell_num):
        Z_count_data = self.sim_model.plasma.abundance[shell_num-1]
        Z_count_data.index.name = 'Z'
        return pd.DataFrame({
            'Element': Z_count_data.index.map(atomic_number2element_symbol),
            # Format to string to show in scientific notation
            'Frac. Ab. (Shell {})'.format(shell_num): Z_count_data.map('{:.6e}'.format)
        })

    def ion_count(self, Z, shell_num):
        ion_num_density = self.sim_model.plasma.ion_number_density[shell_num-1].loc[Z]
        Z_num_density = self.sim_model.plasma.number_density.loc[Z, shell_num-1]
        ion_count_data = ion_num_density/Z_num_density  # Normalization
        ion_count_data.index.name = 'Ion'
        return pd.DataFrame({
            'Species': ion_count_data.index.map(lambda x: species_tuple_to_string((Z, x))),
            'Frac. Ab. (Z={})'.format(Z): ion_count_data.map('{:.6e}'.format, na_action='ignore')
        })

    def level_count(self, ion, Z, shell_num):
        level_num_density = self.sim_model.plasma.level_number_density[shell_num-1].loc[Z, ion]
        ion_num_density = self.sim_model.plasma.ion_number_density[shell_num-1].loc[Z, ion]
        level_count_data = level_num_density/ion_num_density  # Normalization
        level_count_data.index.name = 'Level'
        level_count_data.name = 'Frac. Ab. (Ion {})'.format(ion)
        return level_count_data.map('{:.6e}'.format, na_action='ignore').to_frame()

    
# TODO: Class for model parameters and other stuff

# TODO class for main tab widget