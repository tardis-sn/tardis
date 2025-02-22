from tardis.base import run_tardis
from tardis.io.atom_data.atom_web_download import download_atom_data
from tardis.util.base import (
    atomic_number2element_symbol,
    is_notebook,
    species_tuple_to_string,
)

from tardis.visualization.widgets.util import create_table_widget_shell_info

import pandas as pd
import numpy as np
import panel as pn

pn.extension('tabulator')

class BaseShellInfo:
    """The simulation information that is used by shell info widget"""

    def __init__(
        self,
        t_radiative,
        dilution_factor,
        abundance,
        number_density,
        ion_number_density,
        level_number_density,
    ):
        """Initialize the object with all simulation properties in use

        Parameters
        ----------
        t_radiative : array_like
            Radiative Temperature of each shell of simulation
        dilution_factor : array_like
            Dilution Factor (W) of each shell of simulation model
        abundance : pandas.DataFrame
            Fractional abundance of elements where row labels are atomic number
            and column labels are shell number
        number_density : pandas.DataFrame
            Number densities of elements where row labels are atomic number and
            column labels are shell numbers
        ion_number_density : pandas.DataFrame
            Number densities of ions where rows are multi-indexed with (atomic
            number, ion number) and column labels are shell number
        level_number_density : pandas.DataFrame
            Number densities of levels where rows are multi-indexed with (atomic
            number, ion number, level number) and column labels are shell number
        """
        self.t_radiative = t_radiative
        self.dilution_factor = dilution_factor
        self.abundance = abundance
        self.number_density = number_density
        self.ion_number_density = ion_number_density
        self.level_number_density = level_number_density

    def shells_data(self):
        """Generates shells data in a form that can be used by a table widget

        Returns
        -------
        pandas.DataFrame
            Dataframe containing Rad. Temp. and W against each shell of
            simulation model
        """
        shells_temp_w = pd.DataFrame(
            {
                "Rad. Temp.": self.t_radiative,
                "Dilution Factor": self.dilution_factor,
            }
        )
        shells_temp_w.index = range(
            1, len(self.t_radiative) + 1
        )  # Overwrite index
        shells_temp_w.index.name = "Shell No."
        # Format to string to make qgrid show values in scientific notations
        return shells_temp_w.map(lambda x: f"{x:.6e}")

    def element_count(self, shell_num):
        """Generates fractional abundance of elements present in a specific
        shell in a form that can be used by a table widget

        Parameters
        ----------
        shell_num : int
            Shell number (note: starts from 1, not 0 which is what simulation
            model use)

        Returns
        -------
        pandas.DataFrame
            Dataframe containing element symbol and fractional abundance in a
            specific shell, against each atomic number
        """
        element_count_data = self.abundance[shell_num - 1].copy()
        element_count_data.index.name = "Z"
        element_count_data = element_count_data.fillna(0)
        return pd.DataFrame(
            {
                "Element": element_count_data.index.map(
                    atomic_number2element_symbol
                ),
                # Format to string to show in scientific notation
                f"Frac. Ab. (Shell {shell_num})": element_count_data.map(
                    "{:.6e}".format
                ),
            }
        )

    def ion_count(self, atomic_num, shell_num):
        """Generates fractional abundance of ions of a specific element and
        shell, in a form that can be used by a table widget

        Parameters
        ----------
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number (note: starts from 1, not 0 which is what simulation
            model use)

        Returns
        -------
        pandas.DataFrame
            Dataframe containing ion specie and fractional abundance for a
            specific element, against each ion number
        """
        ion_num_density = self.ion_number_density[shell_num - 1].loc[atomic_num]
        element_num_density = self.number_density.loc[atomic_num, shell_num - 1]
        ion_count_data = ion_num_density / element_num_density  # Normalization
        ion_count_data.index.name = "Ion"
        ion_count_data = ion_count_data.fillna(0)
        return pd.DataFrame(
            {
                "Species": ion_count_data.index.map(
                    lambda x: species_tuple_to_string((atomic_num, x))
                ),
                f"Frac. Ab. (Z={atomic_num})": ion_count_data.map(
                    "{:.6e}".format
                ),
            }
        )

    def level_count(self, ion, atomic_num, shell_num):
        """Generates fractional abundance of levels of a specific ion, element
        and shell, in a form that can be used by a table widget

        Parameters
        ----------
        ion : int
            Ion number (note: starts from 0, same what is used by simulation
            model)
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number (note: starts from 1, not 0 which is what simulation
            model use)

        Returns
        -------
        pandas.DataFrame
            Dataframe containing fractional abundance for a specific ion,
            against each level number
        """
        level_num_density = self.level_number_density[shell_num - 1].loc[
            atomic_num, ion
        ]
        ion_num_density = self.ion_number_density[shell_num - 1].loc[
            atomic_num, ion
        ]
        level_count_data = level_num_density / ion_num_density  # Normalization
        level_count_data.index.name = "Level"
        level_count_data.name = f"Frac. Ab. (Ion={ion})"
        level_count_data = level_count_data.fillna(0)
        return level_count_data.map("{:.6e}".format).to_frame()


class SimulationShellInfo(BaseShellInfo):
    """The simulation information that is used by shell info widget, obtained
    from a TARDIS Simulation object
    """

    def __init__(self, sim_model):
        """Initialize the object with TARDIS Simulation object

        Parameters
        ----------
        sim_model : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
        """
        super().__init__(
            sim_model.simulation_state.t_radiative,
            sim_model.simulation_state.dilution_factor,
            sim_model.simulation_state.abundance,
            sim_model.plasma.number_density,
            sim_model.plasma.ion_number_density,
            sim_model.plasma.level_number_density,
        )


class HDFShellInfo(BaseShellInfo):
    """The simulation information that is used by shell info widget, obtained
    from a simulation HDF file
    """

    def __init__(self, hdf_fpath):
        """Initialize the object with a simulation HDF file

        Parameters
        ----------
        hdf_fpath : str
            A valid path to a simulation HDF file (HDF file must be created
            from a TARDIS Simulation object using :code:`to_hdf` method with
            default arguments)
        """
        with pd.HDFStore(hdf_fpath, "r") as sim_data:
            super().__init__(
                sim_data["/simulation/simulation_state/t_radiative"],
                sim_data["/simulation/simulation_state/dilution_factor"],
                sim_data["/simulation/simulation_state/abundance"],
                sim_data["/simulation/plasma/number_density"],
                sim_data["/simulation/plasma/ion_number_density"],
                sim_data["/simulation/plasma/level_number_density"],
            )


class ShellInfoWidget:
    """The Shell Info Widget to explore abundances in different shells.

    It consists of four interlinked table widgets - shells table; element count,
    ion count and level count tables - allowing to explore fractional abundances
    all the way from elements, to ions, to levels by clicking on the rows of
    tables.
    """

    def __init__(self, shell_info_data):
        """Initialize the object with the shell information of a simulation
        model

        Parameters
        ----------
        shell_info_data : subclass of BaseShellInfo
            Shell information object constructed from Simulation object or HDF
            file
        """
        self.data = shell_info_data

        # Initialize tables with Panel Tabulator
        self.shells_table = create_table_widget_shell_info(
            self.data.shells_data(), [30, 35, 35]
        )

        # Creating the element count table widget
        self.element_count_table = create_table_widget_shell_info(
            self.data.element_count(self.shells_table.value.index[0]),
            [15, 30, 55],
            changeable_col={
                "index": -1,  # since last column will change names
                # Shells table index will give all possible shell numbers
                "other_names": [
                    f"Frac. Ab. (Shell {shell_num})"
                    for shell_num in self.shells_table.value.index
                ],
            },
        )

        # Creating the ion count table widget
        self.ion_count_table = create_table_widget_shell_info(
            self.data.ion_count(
                self.element_count_table.value.index[0],
                self.shells_table.value.index[0],
            ),
            [20, 30, 50],
            changeable_col={
                "index": -1,
                # Since element are same for each shell thus previous table
                # (element counts for shell 1) will give all possible elements
                "other_names": [
                    f"Frac. Ab. (Z={atomic_num})"
                    for atomic_num in self.element_count_table.value.index
                ],
            },
        )

        # Creating the level count table widget
        self.level_count_table = create_table_widget_shell_info(
            self.data.level_count(
                self.ion_count_table.value.index[0],
                self.element_count_table.value.index[0],
                self.shells_table.value.index[0],
            ),
            [30, 70],
            changeable_col={
                "index": -1,
                # Ion values range from 0 to max atomic_num present in
                # element count table
                "other_names": [
                    f"Frac. Ab. (Ion={ion})"
                    for ion in range(
                        0, self.element_count_table.value.index.max() + 1
                    )
                ],
            },
        )

    def update_element_count_table(self, event):
        """Event listener to update the data in the element count table widget 
        based on interaction (row selection event) in the shells table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about the event (see Notes section).

        Notes
        -----
        This function is triggered automatically as an event handler for row 
        selection in the shells table widget. It updates the element count table 
        based on the selected shell number and ensures that the first row in the 
        table is selected to trigger further updates.
        """
        if not self.shells_table.selection:
            return
        shell_num = self.shells_table.selection[0] + 1  # Convert 0-based index to 1-based shell number
        self.element_count_table.value = self.data.element_count(shell_num)
        atomic_num0 = self.element_count_table.value.index[0]
        if self.element_count_table.selection == [0]:
            self.element_count_table.selection = []  # Clear selection to trigger update
        self.element_count_table.selection = [0]  # Select first row

    def update_ion_count_table(self, event):
        """Event listener to update the data in the ion count table widget 
        based on interaction (row selection event) in the element count table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about the event (see Notes section).

        Notes
        -----
        This function is triggered automatically as an event handler for row 
        selection in the element count table widget. It updates the ion count table 
        based on the selected atomic number and shell number and ensures that the 
        first row in the table is selected to trigger further updates.
        """
        if not self.element_count_table.selection or not self.shells_table.selection:
            return
        shell_num = self.shells_table.selection[0] + 1
        atomic_num = self.element_count_table.value.index[self.element_count_table.selection[0]]
        self.ion_count_table.value = self.data.ion_count(atomic_num, shell_num)
        ion0 = self.ion_count_table.value.index[0]
        if self.ion_count_table.selection == [0]:
            self.ion_count_table.selection = []  # Clear selection to trigger update
        self.ion_count_table.selection = [0]  # Select first row

    def update_level_count_table(self, event):
        """Event listener to update the data in the level count table widget 
        based on interaction (row selection event) in the ion count table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about the event (see Notes section).

        Notes
        -----
        This function is triggered automatically as an event handler for row 
        selection in the ion count table widget. It updates the level count table 
        based on the selected ion, atomic number, and shell number.
        """
        if not self.ion_count_table.selection or not self.element_count_table.selection or not self.shells_table.selection:
            return
        shell_num = self.shells_table.selection[0] + 1
        atomic_num = self.element_count_table.value.index[self.element_count_table.selection[0]]
        ion = self.ion_count_table.value.index[self.ion_count_table.selection[0]]
        self.level_count_table.value = self.data.level_count(ion, atomic_num, shell_num)

    def display(
        self,
        shells_table_width=300,
        element_count_table_width=240,
        ion_count_table_width=240,
        level_count_table_width=180,
        **layout_kwargs,
    ):
        """Display the shell info widget by arranging all component widgets 
        and enabling interactions between the table widgets.

        Parameters
        ----------
        shells_table_width : int, optional
            Width of the shells table in pixels, by default 300.
        element_count_table_width : int, optional
            Width of the element count table in pixels, by default 240.
        ion_count_table_width : int, optional
            Width of the ion count table in pixels, by default 240.
        level_count_table_width : int, optional
            Width of the level count table in pixels, by default 180.

        Other Parameters
        ----------------
        **layout_kwargs
            Additional CSS properties to be applied to the table container layout.

        Returns
        -------
        panel.layout.Column
            Shell info widget containing all component tables and descriptive text.
        """
        # Set table widths
        self.shells_table.width = shells_table_width
        self.element_count_table.width = element_count_table_width
        self.ion_count_table.width = ion_count_table_width
        self.level_count_table.width = level_count_table_width

        # Bind event handlers using param.watch
        self.shells_table.param.watch(self.update_element_count_table, 'selection')
        self.element_count_table.param.watch(self.update_ion_count_table, 'selection')
        self.ion_count_table.param.watch(self.update_level_count_table, 'selection')

        # Initial selection
        self.shells_table.selection = [1]  # Start with shell 2 selected (index 1)

        # Layout
        tables_container = pn.Row(
            self.shells_table,
            self.element_count_table,
            self.ion_count_table,
            self.level_count_table,
            styles={'display': 'flex', 'justify-content': 'space-between', **layout_kwargs}
        )
        text = pn.pane.HTML(
            "<b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (i.e all "
            "values sum to 1)<br><b>W</b> denotes <i>Dilution Factor</i> and "
            "<b>Rad. Temp.</b> is <i>Radiative Temperature (in K)</i>"
        )
        return pn.Column(text, tables_container)

def shell_info_from_simulation(sim_model):
    """Create shell info widget from a TARDIS simulation object

    Parameters
    ----------
    sim_model : tardis.simulation.Simulation
        TARDIS Simulation object produced by running a simulation

    Returns
    -------
    ShellInfoWidget
    """
    shell_info_data = SimulationShellInfo(sim_model)
    return ShellInfoWidget(shell_info_data)


def shell_info_from_hdf(hdf_fpath):
    """Create shell info widget from a simulation HDF file

    Parameters
    ----------
    hdf_fpath : str
        A valid path to a simulation HDF file (HDF file must be created
        from a TARDIS Simulation object using :code:`to_hdf` method with
        default arguments)

    Returns
    -------
    ShellInfoWidget
    """
    shell_info_data = HDFShellInfo(hdf_fpath)
    return ShellInfoWidget(shell_info_data)