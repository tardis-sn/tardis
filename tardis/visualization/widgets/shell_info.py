from tardis import run_tardis
from tardis.io.atom_data import download_atom_data
from tardis.util.base import (
    atomic_number2element_symbol,
    species_tuple_to_string,
)

import pandas as pd
import numpy as np
import panel as pn

pn.extension('tabulator')

TABLE_STYLES = {
    'border': '1px solid #ddd',
    'font-size': '12px',
    'background-color': '#fff'
}

TITLE_STYLES = {
    'font-size': '14px',
    'font-weight': 'bold',
    'color': '#333'
}

CONTAINER_STYLES = {
    'padding': '10px',
    'background-color': '#f9f9f9',
    'border': '1px solid #ddd',
    'border-radius': '3px'
}

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
        shells_temp_w.index=range(1, len(self.t_radiative) + 1)

        shells_temp_w.index.name = "Shell No."
        
        return shells_temp_w.map(lambda x: f"{x:.6e}")

    def element_count(self, shell_num):
        """Generates fractional abundance of elements present in a specific shell

        Parameters
        ----------
        shell_num : int
            Shell number starts from 1

        Returns
        -------
        pandas.DataFrame
            Dataframe containing element symbol and fractional abundance
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
        shell

        Parameters
        ----------
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number starts from 1

        Returns
        -------
        pandas.DataFrame
            Dataframe containing ion species and fractional abundance
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
        and shell

        Parameters
        ----------
        ion : int
            Ion number starts from 0
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number starts from 1

        Returns
        -------
        pandas.DataFrame
            Dataframe containing fractional abundance per level,
            
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
    """Shell info
    from a TARDIS Simulation object
    """

    def __init__(self, sim_model):
        """Initialize with a TARDIS Simulation object

        Parameters
        ----------
        sim_model : tardis.simulation.Simulation
            TARDIS Simulation object 
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
    """Shell info from a simulation HDF file
    """

    def __init__(self, hdf_fpath):
        """Initialize with a simulation HDF file

        Parameters
        ----------
        hdf_fpath : str
            Path to a simulation HDF file
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

    def _apply_tabulator_styles(self, tabulator):
        """Apply consistent styles and alignment to a Tabulator widget"""
        tabulator.styles = TABLE_STYLES
        tabulator.text_align = "center"
        tabulator.header_align = "center"
    
         # helper method to create Tabulator widgets with defaults
    def _create_tabulator(self, df, widths, titles=None, **kwargs):
        """Create a Tabulator widget with pre-set common arguements"""
        # defaults = {
        #     "layout": "fit_columns",
        #     "selectable": 1,
        #     "styles": TABLE_STYLES,
        #     "text_align": "center",
        #     "header_align": "center",
        #     # "autosize_mode": "fit_view",
        #     "width": 100,
        # }
        defaults = {
            "layout": "fit_data",
            "selectable": 1,
            "styles": TABLE_STYLES,
            "text_align": "center",
            "header_align": "center",
        }
        if titles:
            defaults["titles"] = titles
        defaults.update(kwargs)
        return pn.widgets.Tabulator(df, widths=widths, **defaults)
    def __init__(self, shell_info_data):
        self.data = shell_info_data

        # Shells table
        shells_df = self.data.shells_data()
        self.shells_title = pn.pane.Markdown("### Shells", margin=(0, 0, 5, 0), styles=TITLE_STYLES)
        self.shells_table = self._create_tabulator(
            shells_df,
            widths={"Shell No.": 80, "Rad. Temp.": 120, "Dilution Factor": 120},
            titles={"Shell No.": "Shell No.", "Rad. Temp.": "Rad. Temp. (K)", "Dilution Factor": "W"}
        )
        
        # Element count table (Shell 1 initially)
        element_df = self.data.element_count(1)
        self.element_title = pn.pane.Markdown("### Elements (Shell 1)", margin=(0, 0, 5, 0), styles=TITLE_STYLES)
        self.element_count_table = self._create_tabulator(
            element_df,
            widths={"Element": 100, "Frac. Ab. (Shell 1)": 140}
        )
        # Ion count table (initially empty, but ensure data is populated on selection)
        self.ion_title = pn.pane.Markdown("### Ions (No Selection)", margin=(0, 0, 5, 0), styles = TITLE_STYLES)
        self.ion_count_table = self._create_tabulator(
            pd.DataFrame(columns=["Species", "Frac. Ab."]),
            widths={"Species": 140, "Frac. Ab.": 140},
        )

        # Level count table (initially empty, but ensure data is populated on selection)
        self.level_title = pn.pane.Markdown("### Levels (No Selection)",
            margin=(0, 0, 5, 0), styles=TITLE_STYLES)
        self.level_count_table = self._create_tabulator(
            pd.DataFrame(columns=["Level", "Frac. Ab."]),
            widths={"Level": 100, "Frac. Ab.": 180},
        )

        # Bind events
        self.shells_table.param.watch(self.update_element_count_table, "selection")
        self.element_count_table.param.watch(self.update_ion_count_table, "selection")
        self.ion_count_table.param.watch(self.update_level_count_table, "selection")

        # Notes with compact styling
        self.notes = pn.pane.HTML(
            """
            <div style="font-size: 12px; padding: 5px; background-color: #f5f5f5; border: 1px solid #ddd; border-radius: 3px; color: #333;">
                <b>Frac. Ab.</b> denotes <i>Fractional Abundances</i> (sum to 1)<br>
                <b>W</b> denotes <i>Dilution Factor</i><br>
                <b>Radiative Temp.</b> is in Kelvin
            </div>
            """,
            width=350,
        )
        self.layout = pn.Column(
            pn.pane.Markdown("# TARDIS Shell Info Explorer", styles={'font-size': '20px', 'font-weight': 'bold', 'margin-bottom': '10px', 'margin-top': '10px', 'color': '#333', 'width': '100%'}),
        self.notes,
            pn.Row(
                pn.Column(self.shells_title, self.shells_table, styles={'padding': '5px', 'border': '1px solid #ff0000'}, sizing_mode="stretch_width"),
                pn.Column(self.element_title, self.element_count_table, styles={'padding': '5px', 'border': '1px solid #9933ff'}, sizing_mode="stretch_width"),
                pn.Column(self.ion_title, self.ion_count_table, styles={'padding': '5px', 'border': '1px solid #0000ff'}, sizing_mode="stretch_width"),
                pn.Column(self.level_title, self.level_count_table, styles={'padding': '5px', 'border': '1px solid #00ff00'}, sizing_mode="stretch_width"),
                styles={'margin': '10px', 'padding': '10px', 'background-color': '#fff', 'border': '1px solid #ddd', 'border-radius': '3px'},
                sizing_mode="scale_both",  # Scales to fit the screen without scrolling
            ),
            sizing_mode="scale_both",
            styles=CONTAINER_STYLES,    
        )
        
        # Initial selection for shells table only if data exists
        if not shells_df.empty:
            self.shells_table.selection = [0]

    def update_element_count_table(self, event):
        """Update element table based on shell selection"""
        selected_rows = event.new
        if selected_rows:
            shell_num = selected_rows[0] + 1
            element_df = self.data.element_count(shell_num)
            self.element_count_table.value = element_df
            self.element_title.object = f"### Elements (Shell {shell_num})"
            self.element_count_table.widths = {
                "Element": 100,
                f"Frac. Ab. (Shell {shell_num})": 140,
            }
            self._apply_tabulator_styles(self.element_count_table)

            if not element_df.empty:
                    self.element_count_table.selection = [0]
            else:
                    self.element_count_table.selection = []
                    self.update_ion_count_table(None)
        
        else:
            self.element_count_table.value = pd.DataFrame(columns=["Element", "Frac. Ab."])
            self.element_title.object = "### Elements (No Shell Selected)"
            # self.element_count_table.selection = []
            self.ion_count_table.value = pd.DataFrame(columns=["Species", "Frac. Ab."])
            self.ion_title.object = "### Ions (No Selection)"
    
            self.level_count_table.value = pd.DataFrame(columns=["Level", "Frac. Ab."])
            self.level_title.object = "### Levels (No Selection)"
            self.update_ion_count_table(None)

            
    def update_ion_count_table(self, event):
        """Update ion table based on element selection"""
        if event and self.shells_table.selection:
            selected_rows = event.new
            shell_num = self.shells_table.selection[0] + 1
            if selected_rows:
                atomic_num = self.element_count_table.value.index[selected_rows[0]]
                ion_df = self.data.ion_count(atomic_num, shell_num)
                self.ion_count_table.value = ion_df
                self.ion_title.object = f"### Ions (Z={atomic_num}, Shell {shell_num})"
                self.ion_count_table.widths = {
                    "Species": 140,
                    f"Frac. Ab. (Z={atomic_num})": 140,
                }
                self._apply_tabulator_styles(self.ion_count_table)

                if not ion_df.empty:
                    self.ion_count_table.selection = [0]
                else:
                    self.ion_count_table.selection = []
                    self.update_level_count_table(None)
            else:
                self.ion_count_table.value = pd.DataFrame(columns=["Species", "Frac. Ab."])
                self.ion_title.object = "### Ions (No Element Selected)"
                self.ion_count_table.selection = []
                self.update_level_count_table(None)
        else:
            self.ion_count_table.value = pd.DataFrame(columns=["Species", "Frac. Ab."])
            self.ion_title.object = "### Ions (No Selection)"
            self.ion_count_table.selection = []
            self.update_level_count_table(None)

    def update_level_count_table(self, event):
        """Update level table based on ion selection"""
        if event and self.shells_table.selection and self.element_count_table.selection:
            selected_rows = event.new
            shell_num = self.shells_table.selection[0] + 1
            atomic_num = self.element_count_table.value.index[self.element_count_table.selection[0]]
            if selected_rows:
                ion = self.ion_count_table.value.index[selected_rows[0]]
                level_df = self.data.level_count(ion, atomic_num, shell_num)
                self.level_count_table.value = level_df
                self.level_title.object = f"### Levels (Ion={ion}, Z={atomic_num}, Shell {shell_num})"
                self.level_count_table.widths = {
                    "Level": 100,
                    f"Frac. Ab. (Ion={ion})": 180,
                }
                self._apply_tabulator_styles(self.level_count_table)

                if not level_df.empty:
                    self.level_count_table.selection = [0]
                else:
                    self.level_count_table.selection = []
            else:
                self.level_count_table.value = pd.DataFrame(columns=["Level", "Frac. Ab."])
                self.level_title.object = "### Levels (No Ion Selected)"
                self.level_count_table.selection = []
        else:
            self.level_count_table.value = pd.DataFrame(columns=["Level", "Frac. Ab."])
            self.level_title.object = "### Levels (No Selection)"
            self.level_count_table.selection = []

    def get_panel(self):
        """Return the Panel object for display"""
        return self.layout


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
