from tardis import run_tardis
from tardis.io.atom_data import download_atom_data
from tardis.util.base import (
    atomic_number2element_symbol,
    species_tuple_to_string,
)
from tardis.util.base import is_notebook
import pandas as pd
import numpy as np
import panel as pn

pn.extension('tabulator')


# Centralized style dictionaries
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
        
    def _normalize_data(self, data):
        """Helper method to normalize data to ensure sum is 1.0
        
        Parameters
        ----------
        data : pandas.Series
            Data to normalize
            
        Returns
        -------
        pandas.Series
            Normalized data that sums to 1.0
        """
        data = data.astype(np.float64)
        if data.sum() > 0:
            data = data / data.sum()
        return data.fillna(0).round(7)

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
        shells_temp_w.index = range(1, len(self.t_radiative) + 1)  # Shells start at 1
        shells_temp_w.index.name = "Shell No."
        return shells_temp_w.map(lambda x: f"{x:.6e}")

    def element_count(self, shell_num, format_for_display=True):
        """Generates fractional abundance of elements present in a specific shell

        Parameters
        ----------
        shell_num : int
            Shell number (starts from 1)

        Returns
        -------
        pandas.DataFrame
            Dataframe with element symbol and fractional abundance
        """
        element_count_data = self.abundance[shell_num - 1].copy()
        element_count_data.index.name = "Z"
        element_count_data = element_count_data.fillna(0)
        df = pd.DataFrame(
            {
                "Element": element_count_data.index.map(atomic_number2element_symbol),
                f"Frac. Ab. (Shell {shell_num})": element_count_data,
            }
        )
        if format_for_display:
            df[f"Frac. Ab. (Shell {shell_num})"] = df[f"Frac. Ab. (Shell {shell_num})"].map("{:.6e}".format)
        return df

    def ion_count(self, atomic_num, shell_num, format_for_display=True):
        """Generates fractional abundance of ions for a specific element and shell

        Parameters
        ----------
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number (starts from 1)

        Returns
        -------
        pandas.DataFrame
            Dataframe with ion species and fractional abundance
        """
        ion_num_density = self.ion_number_density[shell_num - 1].loc[atomic_num]
        element_num_density = self.number_density.loc[atomic_num, shell_num - 1]
        
        # Compute the ratio directly using numpy for consistent results
        ion_count_data = (ion_num_density / element_num_density).astype(np.float64)
        
        # Normalize to ensure sum is 1.0 if values are present
        ion_count_data = self._normalize_data(ion_count_data)
            
        ion_count_data.index.name = "Ion"
        
        df = pd.DataFrame(
            {
                "Species": ion_count_data.index.map(
                    lambda x: species_tuple_to_string((atomic_num, x))
                ),
                f"Frac. Ab. (Z={atomic_num})": ion_count_data,
            }
        )
        if format_for_display:
            df[f"Frac. Ab. (Z={atomic_num})"] = df[f"Frac. Ab. (Z={atomic_num})"].map("{:.6e}".format)
        return df

    def level_count(self, ion, atomic_num, shell_num, format_for_display=True):
        """Generates fractional abundance of levels for a specific ion, element, and shell

        Parameters
        ----------
        ion : int
            Ion number (starts from 0)
        atomic_num : int
            Atomic number of element
        shell_num : int
            Shell number (starts from 1)

        Returns
        -------
        pandas.DataFrame
            Dataframe with fractional abundance per level
        """
        level_num_density = self.level_number_density[shell_num - 1].loc[atomic_num, ion]
        ion_num_density = self.ion_number_density[shell_num - 1].loc[atomic_num, ion]
        
        # Compute the ratio directly using numpy for consistent results
        level_count_data = (level_num_density / ion_num_density).astype(np.float64)
        
        # Normalize to ensure sum is 1.0 if values are present
        level_count_data = self._normalize_data(level_count_data)
            
        level_count_data.index.name = "Level"
        level_count_data.name = f"Frac. Ab. (Ion={ion})"
        
        df = level_count_data.to_frame()
        if format_for_display:
            df[df.columns[0]] = df[df.columns[0]].map("{:.6e}".format)
        return df


class SimulationShellInfo(BaseShellInfo):
    """Shell info from a TARDIS Simulation object"""

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
    """Shell info from a simulation HDF file"""

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
        """
        Apply consistent styles and alignment to a Tabulator widget.

        Parameters
        ----------
        tabulator : panel.widgets.Tabulator
            The tabulator widget to style.
        """
        tabulator.styles = TABLE_STYLES
        tabulator.text_align = "center"
        tabulator.header_align = "center"

    def _apply_styles_to_all_tables(self):
        """
        Apply consistent styles to all tabulator tables in the widget.
        """
        for table in [self.shells_table, self.element_count_table, 
                     self.ion_count_table, self.level_count_table]:
            self._apply_tabulator_styles(table)

    def _reset_table(self, table, title_widget, title_text, columns, selection_callback=None):
        """
        Reset a table to an empty state with default columns.
        
        Parameters
        ----------
        table : panel.widgets.Tabulator
            The table to reset
        title_widget : panel.pane.Markdown
            The title widget to update
        title_text : str
            The new title text
        columns : list
            Column names for the empty dataframe
        selection_callback : function, optional
            Callback to trigger for downstream updates
        """
        table.value = pd.DataFrame(columns=columns)
        title_widget.object = title_text
        table.selection = []
        
        if selection_callback is not None:
            selection_callback(None)

    # helper method to create Tabulator widgets with defaults
    def _create_tabulator(self, df, widths, titles=None, **kwargs):
        """
        Create a Tabulator widget with pre-set common arguments.
        
        This helper method creates a standardized tabulator widget with consistent
        styling and behavior.

        Parameters
        ----------
        df : pandas.DataFrame
            The data to display in the tabulator.
        widths : dict
            Dictionary mapping column names to column widths.
        titles : dict, optional
            Dictionary mapping column names to display titles (default: None).
        **kwargs : dict
            Additional keyword arguments to pass to the Tabulator constructor.

        Returns
        -------
        panel.widgets.Tabulator
            Configured tabulator widget.
        """
        defaults = {
            "layout": "fit_data_table",
            "selectable": 1,
            "styles": TABLE_STYLES,
            "text_align": "center",
            "header_align": "center",
            "height": 600,
        }
        if titles:
            defaults["titles"] = titles
        defaults.update(kwargs)
        return pn.widgets.Tabulator(df, widths=widths, **defaults)

    def _create_table_column(self, title_widget, table_widget, width=260):
        """
        Create a consistent column layout for a table and its title.
        
        Parameters
        ----------
        title_widget : panel.pane.Markdown
            The title widget for the table
        table_widget : panel.widgets.Tabulator
            The table widget
        width : int, optional
            Width of the column in pixels, default: 260
            
        Returns
        -------
        panel.Column
            A column containing the title and table with consistent styling
        """
        return pn.Column(
            title_widget, 
            table_widget, 
            styles={
                'padding': '5px', 
                'border': '1px solid #ddd', 
                'background-color': '#f9f9f9'
            }, 
            width=width
        )

    def __init__(self, shell_info_data):
        """
        Initialize the Shell Info Widget with simulation data.
        
        Sets up the four interconnected tables and establishes event handlers
        for selection changes.

        Parameters
        ----------
        shell_info_data : BaseShellInfo or subclass
            The shell information data source to use for populating tables.
        """
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
            widths={"Element": 100, "Frac. Ab. (Shell 1)": 140},
        )

        # Ion count table (initially empty, but ensure data is populated on selection)
        self.ion_title = pn.pane.Markdown("### Ions (No Selection)", margin=(0, 0, 5, 0), styles=TITLE_STYLES)
        self.ion_count_table = self._create_tabulator(
            pd.DataFrame(columns=["Species", "Frac. Ab."]),
            widths={"Species": 140, "Frac. Ab.": 140},
        )

        # Level count table (initially empty, but ensure data is populated on selection)
        self.level_title = pn.pane.Markdown("### Levels (No Selection)", margin=(0, 0, 5, 0), styles=TITLE_STYLES)
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
            styles={'background-color': '#f5f5f5'}
        )

        # Apply styles to all tables initially
        self._apply_styles_to_all_tables()

        # Improved layout with compact, no-scroll design
        self.layout = pn.Column(
            pn.pane.Markdown(
                "# TARDIS Shell Info Explorer", 
                styles={'font-size': '20px', 'font-weight': 'bold', 'margin-bottom': '10px', 'color': '#333'}
            ),
            self.notes,
            pn.Row(
                self._create_table_column(self.shells_title, self.shells_table),
                self._create_table_column(self.element_title, self.element_count_table),
                self._create_table_column(self.ion_title, self.ion_count_table),
                self._create_table_column(self.level_title, self.level_count_table),
                styles={
                    'margin': '10px', 
                    'padding': '10px', 
                    'background-color': '#fff', 
                    'border': '1px solid #ddd', 
                    'border-radius': '3px'
                },
                sizing_mode="stretch_width", 
            ),
            sizing_mode="stretch_width",
            styles=CONTAINER_STYLES,
        )

        # Initial selection for shells table only if data exists
        if not shells_df.empty:
            self.shells_table.selection = [0]

    def update_element_count_table(self, event):
        """
        Update element table based on shell selection.
        
        This event handler updates the element table when a shell is selected in 
        the shells table. It also triggers downstream updates.

        Parameters
        ----------
        event : param.Event
            The event containing the new selection information.
        """
        selected_rows = event.new if hasattr(event, 'new') else []
        if selected_rows:
            shell_num = selected_rows[0] + 1
            element_df = self.data.element_count(shell_num)

            if isinstance(element_df, pd.DataFrame):
                # Force consistent numeric formatting
                for col in element_df.select_dtypes(include=['float64']).columns:
                    element_df[col] = pd.to_numeric(element_df[col], errors='coerce')
          
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
            self._reset_table(
                self.element_count_table, 
                self.element_title, 
                "### Elements (No Shell Selected)", 
                ["Element", "Frac. Ab."],
                self.update_ion_count_table
            )

    def update_ion_count_table(self, event):
        """
        Update ion table based on element selection.
        
        This event handler updates the ion table when an element is selected in 
        the element count table. It also triggers downstream updates.

        Parameters
        ----------
        event : param.Event or None
            The event containing the new selection information, or None if called directly.
        """
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
                self._reset_table(
                    self.ion_count_table, 
                    self.ion_title, 
                    "### Ions (No Element Selected)", 
                    ["Species", "Frac. Ab."],
                    self.update_level_count_table
                )
        else:
            self._reset_table(
                self.ion_count_table, 
                self.ion_title, 
                "### Ions (No Selection)", 
                ["Species", "Frac. Ab."],
                self.update_level_count_table
            )

    def update_level_count_table(self, event):
        """
        Update level table based on ion selection.
        
        This event handler updates the level table when an ion is selected in 
        the ion count table.

        Parameters
        ----------
        event : param.Event or None
            The event containing the new selection information, or None if called directly.
        """
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
                self._reset_table(
                    self.level_count_table, 
                    self.level_title, 
                    "### Levels (No Ion Selected)", 
                    ["Level", "Frac. Ab."]
                )
        else:
            self._reset_table(
                self.level_count_table, 
                self.level_title, 
                "### Levels (No Selection)", 
                ["Level", "Frac. Ab."]
            )

    def get_panel(self):
        """Return the Panel object for display"""
        return self.layout

    def display(self):
        """
        Display the widget in the current context.
        
        In notebook contexts, returns the Panel object for display.
        In non-notebook contexts, serves the Panel application.

        Returns
        -------
        panel.Column or None
            The main layout if in a notebook context, None otherwise.
        """        
        if is_notebook():
            return self.get_panel()
        else:
            pn.serve(self.get_panel())

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