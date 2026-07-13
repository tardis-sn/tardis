import pandas as pd
import panel as pn
import param

from tardis.util.base import (
    atomic_number2element_symbol,
    species_tuple_to_string,
)
from tardis.util.environment import Environment


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
        # Format to string to show values in scientific notations
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


class ShellInfoWidget(param.Parameterized):
    """The Shell Info Widget to explore abundances in different shells.

    It consists of four interlinked table widgets - shells table; element count,
    ion count and level count tables - allowing to explore fractional abundances
    all the way from elements, to ions, to levels by clicking on the rows of
    tables.
    """

    shell_idx = param.List(default=[], doc="Index of current selected row in shells table")
    atomic_idx = param.List(default=[], doc="Index of current selected row in element count table")
    ion_idx = param.List(default=[], doc="Index of current selected row in ion count table")

    @staticmethod
    def _create_table_widget(data):
        """
        Create a table widget using Panel's Tabulator widget

        Parameters
        ----------
        data: pandas.DataFrame
            Data to be displayed in the table widget

        Returns
        -------
        panel.widgets.Tabulator
        """
        _df = data.copy()

        # Create the table
        table = pn.widgets.Tabulator(
            _df,
            selectable=True,  # Single row selection (radio button style)
            show_index=True,
            sizing_mode='stretch_width',
            height=min(400, max(200, len(_df) * 30 + 50)),
            disabled=True,  # Make cells non-editable
        )
        return table

    def __init__(self, shell_info_data):
        """Initialize the object with the shell information of a simulation
        model

        Parameters
        ----------
        shell_info_data : subclass of BaseShellInfo
            Shell information object constructed from Simulation object or HDF
            file
        """
        super().__init__()
        self.data = shell_info_data

        # Creating the shells data table widget
        self.shells_table = self._create_table_widget(
            self.data.shells_data()
        )

        # Creating the element count table widget
        self.element_count_table = self._create_table_widget(
            self.data.element_count(self.shells_table.value.index[0])
        )

        # Creating the ion count table widget
        self.ion_count_table = self._create_table_widget(
            self.data.ion_count(
                self.element_count_table.value.index[0],
                self.shells_table.value.index[0],
            )
        )

        # Creating the level count table widget
        self.level_count_table = self._create_table_widget(
            self.data.level_count(
                self.ion_count_table.value.index[0],
                self.element_count_table.value.index[0],
                self.shells_table.value.index[0],
            )
        )

        # The indexes will update when user clicks on the rows of tables
        self.shells_table.link(self, selection="shell_idx")
        self.element_count_table.link(self, selection="atomic_idx")
        self.ion_count_table.link(self, selection="ion_idx")

        # Initialize tables.
        self.shells_table.selection = [0]
        self.element_count_table.selection = [0]
        self.ion_count_table.selection = [0]

    @param.depends("shell_idx", watch=True)
    def update_element_count_table(self):
        """Event listener to update the data in element count table widget based
        on interaction (row selected event) in shells table widget.
        """
        # Update element count table data based on selected shell number (shell_idx)
        self.element_count_table.value = self.data.element_count(self.shell_idx[0] + 1)

        self.element_count_table.selection = [0]  # Reset ion_num selection when shell_idx changes

    @param.depends("atomic_idx", "shell_idx", watch=True)
    def update_ion_count_table(self):
        """Event listener to update the data in ion count table widget based
        on interaction (row selected event) in element count table widget.
        """
        # Get the selected shell number and atomic number based on selected rows in shells table and element count table
        shell_num = self.shell_idx[0] + 1
        atomic_num = self.element_count_table.value.index[self.atomic_idx[0]]

        # Update ion count table data based on selected shell number and atomic number
        self.ion_count_table.value = self.data.ion_count(atomic_num, shell_num)

        self.ion_count_table.selection = [0]  # Reset ion count table selection when atomic_idx changes

    @param.depends("ion_idx", "atomic_idx", "shell_idx", watch=True)
    def update_level_count_table(self):
        """Event listener to update the data in level count table widget based
        on interaction (row selected event) in ion count table widget.
        """
        # Get the selected shell number, atomic number and ion number based on selected rows in their tables
        shell_num = self.shell_idx[0] + 1
        atomic_num = self.element_count_table.value.index[self.atomic_idx[0]]
        ion_num = self.ion_count_table.value.index[self.ion_idx[0]]

        self.level_count_table.value = self.data.level_count(
            ion_num, atomic_num, shell_num
        )

    def display(self):
        """Display the shell info widget by putting all component widgets nicely
        together and allowing interaction between the table widgets

        Returns
        -------
        panel.Column
            Shell info widget containing all component widgets
        """
        if not Environment.allows_widget_display():
            print("Please use a notebook to display the widget")
        else:

            # Create Panel layout for the tables
            shell_info_tables_container = pn.Row(
                self.shells_table,
                self.element_count_table,
                self.ion_count_table,
                self.level_count_table,
                sizing_mode='stretch_width'
            )

            # Notes text explaining how to interpret tables widgets' data
            text = pn.pane.HTML(
                "<b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (i.e all "
                "values sum to 1)<br><b>W</b> denotes <i>Dilution Factor</i> and "
                "<b>Rad. Temp.</b> is <i>Radiative Temperature (in K)</i>"
            )

            # Put text vertically before shell info container
            shell_info_widget = pn.Column(text, shell_info_tables_container)
            return shell_info_widget


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
