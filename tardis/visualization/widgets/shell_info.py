from tardis.base import run_tardis
from tardis.io.atom_data.util import download_atom_data
from tardis.util.base import (
    atomic_number2element_symbol,
    species_tuple_to_string,
)

from tardis.visualization.widgets.util import create_table_widget

import pandas as pd
import numpy as np
import ipywidgets as ipw


class BaseShellInfo:
    """The simulation information that is used by shell info widget"""

    def __init__(
        self,
        t_radiative,
        w,
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
        w : array_like
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
        self.w = w
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
            {"Rad. Temp.": self.t_radiative, "W": self.w}
        )
        shells_temp_w.index = range(
            1, len(self.t_radiative) + 1
        )  # Overwrite index
        shells_temp_w.index.name = "Shell No."
        # Format to string to make qgrid show values in scientific notations
        return shells_temp_w.applymap(lambda x: f"{x:.6e}")

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
        element_count_data.fillna(0, inplace=True)
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
        ion_count_data.fillna(0, inplace=True)
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
        level_count_data.fillna(0, inplace=True)
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
            sim_model.model.t_radiative,
            sim_model.model.w,
            sim_model.plasma.abundance,
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
                sim_data["/simulation/model/t_radiative"],
                sim_data["/simulation/model/w"],
                sim_data["/simulation/plasma/abundance"],
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

        # Creating the shells data table widget
        self.shells_table = create_table_widget(
            self.data.shells_data(), [30, 35, 35]
        )

        # Creating the element count table widget
        self.element_count_table = create_table_widget(
            self.data.element_count(self.shells_table.df.index[0]),
            [15, 30, 55],
            changeable_col={
                "index": -1,  # since last column will change names
                # Shells table index will give all possible shell numbers
                "other_names": [
                    f"Frac. Ab. (Shell {shell_num})"
                    for shell_num in self.shells_table.df.index
                ],
            },
        )

        # Creating the ion count table widget
        self.ion_count_table = create_table_widget(
            self.data.ion_count(
                self.element_count_table.df.index[0],
                self.shells_table.df.index[0],
            ),
            [20, 30, 50],
            changeable_col={
                "index": -1,
                # Since element are same for each shell thus previous table
                # (element counts for shell 1) will give all possible elements
                "other_names": [
                    f"Frac. Ab. (Z={atomic_num})"
                    for atomic_num in self.element_count_table.df.index
                ],
            },
        )

        # Creating the level count table widget
        self.level_count_table = create_table_widget(
            self.data.level_count(
                self.ion_count_table.df.index[0],
                self.element_count_table.df.index[0],
                self.shells_table.df.index[0],
            ),
            [30, 70],
            changeable_col={
                "index": -1,
                # Ion values range from 0 to max atomic_num present in
                # element count table
                "other_names": [
                    f"Frac. Ab. (Ion={ion})"
                    for ion in range(
                        0, self.element_count_table.df.index.max() + 1
                    )
                ],
            },
        )

    def update_element_count_table(self, event, qgrid_widget):
        """Event listener to update the data in element count table widget based
        on interaction (row selected event) in shells table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about event (see Notes section)
        qgrid_widget : qgrid.QgridWidget
            QgridWidget instance that fired the event (see Notes section)

        Notes
        -----
        You will never need to pass any of these arguments explicitly. This is
        the expected signature of the function passed to :code:`handler` argument
        of :code:`on` method of a table widget (qgrid.QgridWidget object) as
        explained in `qrid documentation <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.on>`_.
        """
        # Get shell number from row selected in shells_table
        shell_num = event["new"][0] + 1

        # Update data in element_count_table
        self.element_count_table.df = self.data.element_count(shell_num)

        # Get atomic_num of 0th row of element_count_table
        atomic_num0 = self.element_count_table.df.index[0]

        # Also update next table (ion counts) by triggering its event listener
        # Listener won't trigger if last row selected in element_count_table was also 0th
        if self.element_count_table.get_selected_rows() == [0]:
            self.element_count_table.change_selection([])  # Unselect rows
        # Select 0th row in count table which will trigger update_ion_count_table
        self.element_count_table.change_selection([atomic_num0])

    def update_ion_count_table(self, event, qgrid_widget):
        """Event listener to update the data in ion count table widget based
        on interaction (row selected event) in element count table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about event (see Notes section)
        qgrid_widget : qgrid.QgridWidget
            QgridWidget instance that fired the event (see Notes section)

        Notes
        -----
        You will never need to pass any of these arguments explicitly. This is
        the expected signature of the function passed to :code:`handler` argument
        of :code:`on` method of a table widget (qgrid.QgridWidget object) as
        explained in `qrid documentation <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.on>`_.
        """
        # Don't execute function if no row was selected, implicitly i.e. by api
        if event["new"] == [] and event["source"] == "api":
            return

        # Get shell no. & atomic_num from rows selected in previous tables
        shell_num = self.shells_table.get_selected_rows()[0] + 1
        atomic_num = self.element_count_table.df.index[event["new"][0]]

        # Update data in ion_count_table
        self.ion_count_table.df = self.data.ion_count(atomic_num, shell_num)

        # Also update next table (level counts) by triggering its event listener
        ion0 = self.ion_count_table.df.index[0]
        if self.ion_count_table.get_selected_rows() == [0]:
            self.ion_count_table.change_selection([])
        self.ion_count_table.change_selection([ion0])

    def update_level_count_table(self, event, qgrid_widget):
        """Event listener to update the data in level count table widget based
        on interaction (row selected event) in ion count table widget.

        Parameters
        ----------
        event : dict
            Dictionary that holds information about event (see Notes section)
        qgrid_widget : qgrid.QgridWidget
            QgridWidget instance that fired the event (see Notes section)

        Notes
        -----
        You will never need to pass any of these arguments explicitly. This is
        the expected signature of the function passed to :code:`handler` argument
        of :code:`on` method of a table widget (qgrid.QgridWidget object) as
        explained in `qrid documentation <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.on>`_.
        """
        # Don't execute function if no row was selected implicitly (by api)
        if event["new"] == [] and event["source"] == "api":
            return

        # Get shell no., atomic_num, ion from selected rows in previous tables
        shell_num = self.shells_table.get_selected_rows()[0] + 1
        atomic_num = self.element_count_table.df.index[
            self.element_count_table.get_selected_rows()[0]
        ]
        ion = self.ion_count_table.df.index[event["new"][0]]

        # Update data in level_count_table
        self.level_count_table.df = self.data.level_count(
            ion, atomic_num, shell_num
        )

    def display(
        self,
        shells_table_width="30%",
        element_count_table_width="24%",
        ion_count_table_width="24%",
        level_count_table_width="18%",
        **layout_kwargs,
    ):
        """Display the shell info widget by putting all component widgets nicely
        together and allowing interaction between the table widgets

        Parameters
        ----------
        shells_table_width : str, optional
            CSS :code:`width` property value for shells table, by default '30%'
        element_count_table_width : str, optional
            CSS :code:`width` property value for element count table, by default '24%'
        ion_count_table_width : str, optional
            CSS :code:`width` property value for ion count table, by default '24%'
        level_count_table_width : str, optional
            CSS :code:`width` property value for level count table, by default '18%'

        Other Parameters
        ----------------
        **layout_kwargs
            Any valid CSS properties to be passed to the :code:`layout` attribute
            of table widgets container (HTML :code:`div`) as explained in
            `ipywidgets documentation <https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20Styling.html#The-layout-attribute>`_

        Returns
        -------
        ipywidgets.Box
            Shell info widget containing all component widgets
        """
        # CSS properties of the layout of shell info tables container
        tables_container_layout = dict(
            display="flex",
            align_items="flex-start",
            justify_content="space-between",
        )
        tables_container_layout.update(layout_kwargs)

        # Setting tables' widths
        self.shells_table.layout.width = shells_table_width
        self.element_count_table.layout.width = element_count_table_width
        self.ion_count_table.layout.width = ion_count_table_width
        self.level_count_table.layout.width = level_count_table_width

        # Attach event listeners to table widgets
        self.shells_table.on(
            "selection_changed", self.update_element_count_table
        )
        self.element_count_table.on(
            "selection_changed", self.update_ion_count_table
        )
        self.ion_count_table.on(
            "selection_changed", self.update_level_count_table
        )

        # Putting all table widgets in a container styled with tables_container_layout
        shell_info_tables_container = ipw.Box(
            [
                self.shells_table,
                self.element_count_table,
                self.ion_count_table,
                self.level_count_table,
            ],
            layout=ipw.Layout(**tables_container_layout),
        )
        self.shells_table.change_selection([1])

        # Notes text explaining how to interpret tables widgets' data
        text = ipw.HTML(
            "<b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (i.e all "
            "values sum to 1)<br><b>W</b> denotes <i>Dilution Factor</i> and "
            "<b>Rad. Temp.</b> is <i>Radiative Temperature (in K)</i>"
        )

        # Put text horizontally before shell info container
        shell_info_widget = ipw.VBox([text, shell_info_tables_container])
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
