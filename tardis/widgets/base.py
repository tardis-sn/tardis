from tardis.base import run_tardis
from tardis.io.atom_data.util import download_atom_data
from tardis.util.base import atomic_number2element_symbol, species_tuple_to_string

import pandas as pd
import numpy as np
import qgrid
import ipywidgets as ipw

# TODO: Fetch saved Simulation object from memory
sim = run_tardis('tardis_example.yml')  # TODO: Make the file path work


class ShellInfoData():
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


class ShellInfoWidget():
    def __init__(self):
        # TODO: Chunk-ify the code in more functions

        # Setting the layout options to be used by table widgets
        grid_options = {
            'sortable': False,
            'filterable': False,
            'editable': False,
            'minVisibleRows': 2,
            'maxVisibleRows': max_rows
        }
        column_options = {
            'minWidth': None,
        }
        # Since forceFitColumns is enabled by default in grid_options,
        # the column widths (when all specified) get applied in proportions,
        # despite their original unit is px (all values below sum to 100)
        shells_col_widths = [30, 35, 35]
        Z_count_col_widths = [15, 30, 55]
        ion_count_col_widths = [20, 30, 50]
        level_count_col_widths = [30, 70]

        def column_widths_definitions(df, col_widths):
            '''
            Generate column definition dictionary from the widths specified for 
            each column (col_widths) including index as a column, in a dataframe (df)
            '''
            cols_with_index = [df.index.name] + df.columns.to_list()
            return {col_name: {'width': col_width}
                    for col_name, col_width in zip(cols_with_index, col_widths)}

        # Creating the shells data table widget
        self.shells_table = qgrid.show_grid(shells_data,
                                            grid_options=grid_options,
                                            column_options=column_options,
                                            column_definitions=column_widths_definitions(
                                                shells_data,
                                                shells_col_widths
                                            ))

        # Creating the Z count table widget
        Z_count_shell1 = Z_count(1)
        Z_column_widths_definitions = column_widths_definitions(
            Z_count_shell1,
            Z_count_col_widths
        )
        for shell_num in range(1, 21):
            Z_column_widths_definitions['Frac. Ab. (Shell {})'.format(
                shell_num)] = {'width': Z_count_col_widths[-1]}
        self.Z_count_table = qgrid.show_grid(Z_count_shell1,
                                             grid_options=grid_options,
                                             column_options=column_options,
                                             column_definitions=Z_column_widths_definitions)

        # Creating the ion count table widget
        ion_count_Z8_shell1 = ion_count(8, 1)
        ion_column_widths_definitions = column_widths_definitions(
            ion_count_Z8_shell1,
            ion_count_col_widths
        )
        for Z in Z_count_shell1.index:
            ion_column_widths_definitions['Frac. Ab. (Z={})'.format(
                Z)] = {'width': ion_count_col_widths[-1]}
        self.ion_count_table = qgrid.show_grid(ion_count_Z8_shell1,
                                               grid_options=grid_options,
                                               column_options=column_options,
                                               column_definitions=ion_column_widths_definitions)

        # Creating the level count table widget
        level_count_ion0_Z8_shell1 = level_count(0, 8, 1)
        level_column_widths_definitions = column_widths_definitions(
            level_count_ion0_Z8_shell1,
            level_count_col_widths
        )
        for level in range(0, 21):
            level_column_widths_definitions['Frac. Ab. (Ion {})'.format(
                level)] = {'width': level_count_col_widths[-1]}
        self.level_count_table = qgrid.show_grid(level_count_ion0_Z8_shell1,
                                                 grid_options=grid_options,
                                                 column_options=column_options,
                                                 column_definitions=level_column_widths_definitions)

    def update_Z_count_table(self, event, qgrid_widget):
        # Get shell number from row selected in shells_table
        shell_num = event['new'][0]+1

        # Update data in Z_count_table
        self.Z_count_table.df = Z_count(shell_num)

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
        self.ion_count_table.df = ion_count(Z, shell_num)

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
        self.level_count_table.df = level_count(ion, Z, shell_num)

    def display(tables_container_layout=ipw.Layout(display='flex',
                                                   align_items='flex-start',
                                                   justify_content='space-between'),
                shells_table_width='32%',
                Z_count_table_width='24%',
                ion_count_table_width='24%',
                level_count_table_width='18%',
                model_param_info=modelparam
                ):

        # Setting tables' widths
        shells_table.layout.width = shells_table_width
        Z_count_table.layout.width = Z_count_table_width
        ion_count_table.layout.width = ion_count_table_width
        level_count_table.layout.width = level_count_table_width

        # Attach event listeners to tables
        self.shells_table.on('selection_changed', update_Z_count_table)
        self.Z_count_table.on('selection_changed', update_ion_count_table)
        self.ion_count_table.on('selection_changed', update_level_count_table)

        # Putting all tables in a container styled with tables_container_layout
        shell_info_tables_container = ipw.Box(
            [shells_table, Z_count_table, ion_count_table, level_count_table],
            layout=tables_container_layout)
        shells_table.change_selection([1])

        # Key to Abbreviation text
        text = ipw.HTML(
            '<b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (i.e all values sum to 1)<br>'
            '<b>W</b> denotes <i>Dilution Factor</i> and <b>Rad. Temp.</b> is <i>Radiative Temperatire (in K)</i>'
        )

        # Put text horizontally before shell info container
        shell_info_widget = ipw.VBox([text, shell_info_tables_container])
        return shell_info_widget

# TODO: Class for model parameters and other stuff

# TODO class for main tab widget
