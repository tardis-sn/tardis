from tardis.base import run_tardis
from tardis.io.atom_data.util import download_atom_data
from tardis.util.base import atomic_number2element_symbol, species_tuple_to_string
from tardis.simulation import Simulation

import pandas as pd
import numpy as np
import qgrid
import ipywidgets as ipw


class BaseShellInfo():
    def __init__(self, t_radiative, w, abundance, number_density,
                 ion_number_density, level_number_density):
        self.t_radiative = t_radiative
        self.w = w
        self.abundance = abundance
        self.number_density = number_density
        self.ion_number_density = ion_number_density
        self.level_number_density = level_number_density

    def shells_data(self):
        shells_temp_w = pd.DataFrame({'Rad. Temp.': self.t_radiative,
                                      'W': self.w})
        shells_temp_w.index = range(
            1, len(self.t_radiative)+1)  # Overwrite index
        shells_temp_w.index.name = 'Shell No.'
        # Format to string to make qgrid show values in scientific notations
        return shells_temp_w.applymap(lambda x: '{:.6e}'.format(x))

    def Z_count(self, shell_num):
        Z_count_data = self.abundance[shell_num-1]
        Z_count_data.index.name = 'Z'
        return pd.DataFrame({
            'Element': Z_count_data.index.map(atomic_number2element_symbol),
            # Format to string to show in scientific notation
            'Frac. Ab. (Shell {})'.format(shell_num): Z_count_data.map('{:.6e}'.format)
        })

    def ion_count(self, Z, shell_num):
        ion_num_density = self.ion_number_density[shell_num-1].loc[Z]
        Z_num_density = self.number_density.loc[Z, shell_num-1]
        ion_count_data = ion_num_density/Z_num_density  # Normalization
        ion_count_data.index.name = 'Ion'
        return pd.DataFrame({
            'Species': ion_count_data.index.map(lambda x: species_tuple_to_string((Z, x))),
            'Frac. Ab. (Z={})'.format(Z): ion_count_data.map('{:.6e}'.format, na_action='ignore')
        })

    def level_count(self, ion, Z, shell_num):
        level_num_density = self.level_number_density[shell_num-1].loc[Z, ion]
        ion_num_density = self.ion_number_density[shell_num-1].loc[Z, ion]
        level_count_data = level_num_density/ion_num_density  # Normalization
        level_count_data.index.name = 'Level'
        level_count_data.name = 'Frac. Ab. (Ion={})'.format(ion)
        return level_count_data.map('{:.6e}'.format, na_action='ignore').to_frame()


class SimulationShellInfo(BaseShellInfo):
    def __init__(self, sim_model):
        super().__init__(
            sim_model.model.t_radiative,
            sim_model.model.w,
            sim_model.plasma.abundance,
            sim_model.plasma.number_density,
            sim_model.plasma.ion_number_density,
            sim_model.plasma.level_number_density
        )


class HDFShellInfo(BaseShellInfo):
    def __init__(self, hdf_fpath):
        with pd.HDFStore(hdf_fpath, 'r') as sim_data:
            super().__init__(
                sim_data['/simulation/model/t_radiative'],
                sim_data['/simulation/model/w'],
                sim_data['/simulation/plasma/abundance'],
                sim_data['/simulation/plasma/number_density'],
                sim_data['/simulation/plasma/ion_number_density'],
                sim_data['/simulation/plasma/level_number_density']
            )


class ShellInfoWidget():
    def __init__(self, sim_obj_or_path):
        if isinstance(sim_obj_or_path, Simulation):
            self.data = SimulationShellInfo(sim_obj_or_path)
        elif isinstance(sim_obj_or_path, str):
            self.data = HDFShellInfo(sim_obj_or_path)
        else:
            raise TypeError(
                "Passed argument is of invalid type: {}\nOnly "
                "tardis.simulation.Simulation object or a path to simulation "
                "HDF file (str object) is allowed!".format(type(sim_obj_or_path)))

        # Creating the shells data table widget
        self.shells_table = self.create_table_widget(
            self.data.shells_data(),
            [30, 35, 35]
        )

        # Creating the Z count table widget
        self.Z_count_table = self.create_table_widget(
            self.data.Z_count(self.shells_table.df.index[0]),
            [15, 30, 55],
            -1,  # since last column will change names
            # Shells table index will give all possible shell numbers
            ['Frac. Ab. (Shell {})'.format(shell_num)
             for shell_num in self.shells_table.df.index]
        )

        # Creating the ion count table widget
        self.ion_count_table = self.create_table_widget(
            self.data.ion_count(self.Z_count_table.df.index[0],
                                self.shells_table.df.index[0]),
            [20, 30, 50],
            -1,
            # Since Z are same for each shell thus previous table (Z counts
            # for shell 1) will give all possible Z
            ['Frac. Ab. (Z={})'.format(Z) for Z in self.Z_count_table.df.index]
        )

        # Creating the level count table widget
        self.level_count_table = self.create_table_widget(
            self.data.level_count(self.ion_count_table.df.index[0],
                                  self.Z_count_table.df.index[0],
                                  self.shells_table.df.index[0]),
            [30, 70],
            -1,
            # Ion values range from 0 to maximum Z present in Z counts table
            ['Frac. Ab. (Ion={})'.format(ion)
             for ion in range(0, self.Z_count_table.df.index.max()+1)]
        )

    def create_table_widget(self, data, col_widths, changeable_col_idx=None,
                            other_col_names=None):
        # Setting the options to be used for creating table widgets
        grid_options = {
            'sortable': False,
            'filterable': False,
            'editable': False,
            'minVisibleRows': 2
        }
        column_options = {
            'minWidth': None,
        }

        # Check whether passed col_widths list is correct or not
        if len(col_widths) != data.shape[1]+1:
            raise ValueError('Size of column widths list do not match with '
                             'number of columns + 1 (index) in dataframe')

        # Note: Since forceFitColumns is enabled by default in grid_options,
        # the column widths (when all specified) get applied in proportions,
        # despite their original unit is px thus it's better they sum to 100
        if sum(col_widths) != 100:
            raise ValueError('Column widths are not proportions of 100 (i.e. '
                             'they do not sum to 100)')

        # Preparing dictionary that defines column widths
        cols_with_index = [data.index.name] + data.columns.to_list()
        column_widths_definitions = {col_name: {'width': col_width}
                                     for col_name, col_width in zip(cols_with_index, col_widths)}

        # We also need to define widths for different names of changeable column
        if changeable_col_idx:
            column_widths_definitions.update(
                {col_name: {'width': col_widths[changeable_col_idx]}
                 for col_name in other_col_names})

        # Create the table widget using qgrid
        return qgrid.show_grid(data,
                               grid_options=grid_options,
                               column_options=column_options,
                               column_definitions=column_widths_definitions)

    def update_Z_count_table(self, event, qgrid_widget):
        # Get shell number from row selected in shells_table
        shell_num = event['new'][0]+1

        # Update data in Z_count_table
        self.Z_count_table.df = self.data.Z_count(shell_num)

        # Get Z of 0th row of Z_count_table
        Z0 = self.Z_count_table.df.index[0]

        # Also update next table (ion counts) by triggering its event listener
        # Listener won't trigger if last row selected in Z_count_table was also 0th
        if self.Z_count_table.get_selected_rows() == [0]:
            self.Z_count_table.change_selection([])  # Unselect rows
        # Select 0th row in count table which will trigger update_ion_count_table
        self.Z_count_table.change_selection([Z0])

    def update_ion_count_table(self, event, qgrid_widget):
        # Don't execute function if no row was selected, implicitly i.e. by api
        if event['new'] == [] and event['source'] == 'api':
            return

        # Get shell no. & Z from rows selected in previous tables
        shell_num = self.shells_table.get_selected_rows()[0]+1
        Z = self.Z_count_table.df.index[event['new'][0]]

        # Update data in ion_count_table
        self.ion_count_table.df = self.data.ion_count(Z, shell_num)

        # Also update next table (level counts) by triggering its event listener
        ion0 = self.ion_count_table.df.index[0]
        if self.ion_count_table.get_selected_rows() == [0]:
            self.ion_count_table.change_selection([])
        self.ion_count_table.change_selection([ion0])

    def update_level_count_table(self, event, qgrid_widget):
        # Don't execute function if no row was selected implicitly (by api)
        if event['new'] == [] and event['source'] == 'api':
            return

        # Get shell no., Z, ion from selected rows in previous tables
        shell_num = self.shells_table.get_selected_rows()[0]+1
        Z = self.Z_count_table.df.index[self.Z_count_table.get_selected_rows()[
            0]]
        ion = self.ion_count_table.df.index[event['new'][0]]

        # Update data in level_count_table
        self.level_count_table.df = self.data.level_count(ion, Z, shell_num)

    def display(self,
                tables_container_layout=ipw.Layout(display='flex',
                                                   align_items='flex-start',
                                                   justify_content='space-between'),
                shells_table_width='30%',
                Z_count_table_width='24%',  # 25%
                ion_count_table_width='24%',  # 25%
                level_count_table_width='18%'
                ):

        # Setting tables' widths
        self.shells_table.layout.width = shells_table_width
        self.Z_count_table.layout.width = Z_count_table_width
        self.ion_count_table.layout.width = ion_count_table_width
        self.level_count_table.layout.width = level_count_table_width

        # Attach event listeners to tables
        self.shells_table.on('selection_changed', self.update_Z_count_table)
        self.Z_count_table.on('selection_changed', self.update_ion_count_table)
        self.ion_count_table.on('selection_changed',
                                self.update_level_count_table)

        # Putting all tables in a container styled with tables_container_layout
        shell_info_tables_container = ipw.Box(
            [self.shells_table, self.Z_count_table,
                self.ion_count_table, self.level_count_table],
            layout=tables_container_layout)
        self.shells_table.change_selection([1])

        # Notes text explaining how to interpret tables
        text = ipw.HTML(
            '<b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (i.e all '
            'values sum to 1)<br><b>W</b> denotes <i>Dilution Factor</i> and '
            '<b>Rad. Temp.</b> is <i>Radiative Temperature (in K)</i>'
        )

        # Put text horizontally before shell info container
        shell_info_widget = ipw.VBox([text, shell_info_tables_container])
        return shell_info_widget

# TODO: Class for model parameters and other stuff

# TODO class for main tab widget
