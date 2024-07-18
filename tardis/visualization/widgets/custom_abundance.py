"""Class to create and display Custom Abundance Widget."""
import os
import yaml
import numpy as np
import pandas as pd
import ipywidgets as ipw
import plotly.graph_objects as go
from astropy import units as u
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z
from pathlib import Path

import tardis
from tardis.io.model.readers.generic_readers import read_uniform_mass_fractions
from tardis.util.base import (
    quantity_linspace,
    is_valid_nuclide_or_elem,
    is_notebook,
)
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.io.model.parse_density_configuration import (
    calculate_power_law_density,
    calculate_exponential_density,
)
from tardis.io.atom_data.base import AtomData
from tardis.io.configuration.config_validator import validate_dict
from tardis.io.model.readers.csvy import load_csvy
from tardis.io.model.readers.csvy import (
    parse_csv_mass_fractions,
)
from tardis.util.base import atomic_number2element_symbol, quantity_linspace
from tardis.visualization.tools.convergence_plot import transition_colors
from tardis.visualization.widgets.util import debounce

BASE_DIR = tardis.__path__[0]
YAML_DELIMITER = "---"


class CustomAbundanceWidgetData:
    """The model information and data that required in custom
    abundance widget.

    Attributes
    ----------
    elements : list of str
        A list of elements or isotopes' symbols.
    """

    def __init__(self, density_t_0, density, abundance, velocity):
        """Initialize CustomAbundanceWidgetData with model information.

        Parameters
        ----------
        density_t_0 : astropy.units.quantity.Quantity
            Initial time for the density in the model.
        density : astropy.units.quantity.Quantity
        abundance : pd.DataFrame
        velocity : astropy.units.quantity.Quantity
        """
        self.density_t_0 = density_t_0.to("day")
        self.density = density.to("g cm^-3")
        self.abundance = abundance
        self.velocity = velocity.to("km/s")
        self.elements = self.get_symbols()

    def get_symbols(self):
        """Get symbol string from atomic number and mass number."""
        str_symbols = np.array(
            self.abundance.index.get_level_values(0).map(
                atomic_number2element_symbol
            )
        )
        str_mass = np.array(
            self.abundance.index.get_level_values(1), dtype="str"
        )
        return np.add(str_symbols, str_mass)

    @classmethod
    def from_csvy(cls, fpath):
        """Create a new CustomAbundanceWidgetData instance with data
        from CSVY file.

        Parameters
        ----------
        fpath : str
            the path of CSVY file.

        Returns
        -------
        CustomAbundanceWidgetData
        """
        csvy_model_config, csvy_model_data = load_csvy(fpath)
        csvy_schema_file = os.path.join(
            BASE_DIR, "io", "schemas", "csvy_model.yml"
        )
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_file)
        )

        if hasattr(csvy_model_config, "velocity"):
            velocity = quantity_linspace(
                csvy_model_config.velocity.start,
                csvy_model_config.velocity.stop,
                csvy_model_config.velocity.num + 1,
            ).cgs
        else:
            velocity_field_index = [
                field["name"] for field in csvy_model_config.datatype.fields
            ].index("velocity")
            velocity_unit = u.Unit(
                csvy_model_config.datatype.fields[velocity_field_index]["unit"]
            )
            velocity = csvy_model_data["velocity"].values * velocity_unit

        no_of_shells = len(velocity) - 1

        if hasattr(csvy_model_config, "density"):
            adjusted_velocity = velocity.insert(0, 0)
            v_middle = (
                adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
            )
            no_of_shells = len(adjusted_velocity) - 1

            d_conf = csvy_model_config.density
            density_type = d_conf.type
            if density_type == "branch85_w7":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.w7_v_0, d_conf.w7_rho_0, -7
                )
                time_0 = d_conf.w7_time_0
            elif density_type == "uniform":
                density_0 = d_conf.value.to("g cm^-3") * np.ones(no_of_shells)
                time_0 = d_conf.get("time_0", 0 * u.day)
            elif density_type == "power_law":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.v_0, d_conf.rho_0, d_conf.exponent
                )
                time_0 = d_conf.get("time_0", 0 * u.day)
            elif density_type == "exponential":
                density_0 = calculate_exponential_density(
                    v_middle, d_conf.v_0, d_conf.rho_0
                )
                time_0 = d_conf.get("time_0", 0 * u.day)
            else:
                raise ValueError(f"Unrecognized density type " f"{d_conf.type}")
        else:
            density_field_index = [
                field["name"] for field in csvy_model_config.datatype.fields
            ].index("density")
            density_unit = u.Unit(
                csvy_model_config.datatype.fields[density_field_index]["unit"]
            )
            density_0 = csvy_model_data["density"].values * density_unit
            time_0 = csvy_model_config.model_density_time_0

        if hasattr(csvy_model_config, "abundance"):
            abundances_section = csvy_model_config.abundance
            abundance, isotope_abundance = read_uniform_mass_fractions(
                abundances_section, no_of_shells
            )
        else:
            _, abundance, isotope_abundance = parse_csv_mass_fractions(
                csvy_model_data
            )
            abundance = abundance.loc[:, 1:]
            abundance.columns = np.arange(abundance.shape[1])
            isotope_abundance = isotope_abundance.loc[:, 1:]
            isotope_abundance.columns = np.arange(isotope_abundance.shape[1])

        abundance = abundance.replace(np.nan, 0.0)
        abundance = abundance[abundance.sum(axis=1) > 0]
        isotope_abundance = isotope_abundance.replace(np.nan, 0.0)
        isotope_abundance = isotope_abundance[isotope_abundance.sum(axis=1) > 0]

        # Combine elements and isotopes to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        return cls(
            density_t_0=time_0,
            density=density_0,
            abundance=abundance,
            velocity=velocity,
        )

    @classmethod
    def from_yml(cls, fpath, atom_data=None):
        """Create a new CustomAbundanceWidgetData instance with data
        from YAML file.

        Parameters
        ----------
        fpath : str
            The path of YAML file.

        Returns
        -------
        CustomAbundanceWidgetData
        """
        config = Configuration.from_yaml(fpath)
        if atom_data is None:
            atom_data = AtomData.from_hdf(config.atom_data)
        if hasattr(config, "csvy_model"):
            simulation_state = SimulationState.from_csvy(
                config, atom_data=atom_data
            )
        else:
            simulation_state = SimulationState.from_config(
                config, atom_data=atom_data
            )

        velocity = simulation_state.velocity
        density_t_0 = simulation_state.time_explosion
        density = simulation_state.density
        abundance = simulation_state.abundance
        isotopic_mass_fraction = (
            simulation_state.composition.isotopic_mass_fraction
        )

        # Combine elements and isotopes to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotopic_mass_fraction])
        abundance.sort_index(inplace=True)

        return cls(
            density_t_0=density_t_0,
            density=density,
            abundance=abundance,
            velocity=velocity,
        )

    @classmethod
    def from_hdf(cls, fpath):
        """Create a new CustomAbundanceWidgetData instance with data
        from HDF file.

        Parameters
        ----------
        fpath : str
            the path of HDF file.

        Returns
        -------
        CustomAbundanceWidgetData
        """
        with pd.HDFStore(fpath, "r") as hdf:
            abundance = hdf["/simulation/plasma/abundance"]
            _density_t_0 = hdf["/simulation/model/homologous_density/scalars"]
            _density = hdf["/simulation/model/homologous_density/density_0"]
            v_inner = hdf["/simulation/model/v_inner"]
            v_outer = hdf["/simulation/model/v_outer"]

        density_t_0 = float(_density_t_0) * u.s
        density = np.array(_density) * u.g / (u.cm) ** 3
        velocity = np.append(v_inner, v_outer[len(v_outer) - 1]) * u.cm / u.s

        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)

        return cls(
            density_t_0=density_t_0,
            density=density,
            abundance=abundance,
            velocity=velocity,
        )

    @classmethod
    def from_simulation(cls, sim):
        """Create a new CustomAbundanceWidgetData instance from a
        Simulation object.

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        CustomAbundanceWidgetData
        """
        abundance = sim.simulation_state.abundance.copy()
        isotope_abundance = (
            sim.simulation_state.composition.raw_isotope_abundance.copy()
        )

        # integrate element and isotope to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        velocity = sim.simulation_state.velocity
        density_t_0 = sim.simulation_state.time_explosion
        density = sim.simulation_state.density

        return cls(
            density_t_0=density_t_0,
            density=density,
            abundance=abundance,
            velocity=velocity,
        )


class CustomYAML(yaml.YAMLObject):
    """A custom YAML object generated by required properties."""

    def __init__(
        self, name, d_time_0, i_time_0, v_inner_boundary, v_outer_boundary
    ):
        """Initialize CustomYAML object with given properties.

        Parameters
        ----------
        name : str
            Name of the YAML file.
        d_time_0 : astropy.units.quantity.Quantity
            Initial time for the density in the model.
        i_time_0 : astropy.units.quantity.Quantity
            Initial time for isotope decay. Set to 0 for no isotopes.
        v_inner_boundary : astropy.units.quantity.Quantity
            Velocity of the inner boundary.
        v_outer_boundary : astropy.units.quantity.Quantity
            Velocity of the outer boundary.
        """
        self.name = name
        self.model_density_time_0 = d_time_0
        self.model_isotope_time_0 = i_time_0
        self.tardis_model_config_version = "v1.0"
        self.datatype = {}
        self.datatype["fields"] = []
        self.v_inner_boundary = v_inner_boundary
        self.v_outer_boundary = v_outer_boundary

    def create_fields_dict(self, elements):
        """Create a dictionary to store the items in 'fields' part.

        Parameters
        ----------
        elements : list of str
            A list of elements or isotopes' symbols.
        """
        for i in range(len(elements) + 2):
            field = {}

            if i == 0:
                field["name"] = "velocity"
                field["unit"] = "km/s"
            elif i == 1:
                field["name"] = "density"
                field["unit"] = "g/cm^3"
            else:
                field["name"] = elements[i - 2]
                field["desc"] = f"fractional {elements[i-2]} abundance"

            self.datatype["fields"].append(field)


class CustomAbundanceWidget:
    """Widget to edit abundances and densities of simulation model
    graphically.

    It generates a GUI based on input data. The GUI has a plot section
    to visualize the profile, an edit section to allow the user directly
    edit abundance and density profile, and an output section to output
    the model to CSVY file.

    Attributes
    ----------
    shell_no : int
        The selected shell number.
    no_of_shells : int
        The number of shells in the model.
    no_of_elements : int
        The number of elements in the model.
    checked_list : list of bool
        A list of bool to record whether the checkbox is checked.
        The index of the bool corresponds to the index of checkbox.
    fig : plotly.graph_objs._figurewidget.FigureWidget
        The figure object of abundance density plot.
    plot_cmap : str, default: "jet", optional
        String defines the colormap used in abundance density plot.
    _trigger : bool
        If False, disable the callback when abundance input is changed.
    """

    error_view = ipw.Output()

    def __init__(self, widget_data):
        """Initialize CustomAbundanceWidget with data and generate
        the widgets and plot.

        Parameters
        ----------
        widget_data : CustomAbundanceWidgetData
        """
        self.data = widget_data
        self.fig = go.FigureWidget()
        self._trigger = True

        self.create_widgets()
        self.generate_abundance_density_plot()
        self.density_editor = DensityEditor(
            self.data,
            self.fig,
            self.dpd_shell_no,
        )

    @property
    def shell_no(self):
        return self.dpd_shell_no.value

    @shell_no.setter
    def shell_no(self, value):
        self.dpd_shell_no.value = value

    @property
    def no_of_shells(self):
        return self.data.abundance.shape[1]

    @property
    def no_of_elements(self):
        return self.data.abundance.shape[0]

    @property
    def checked_list(self):  # A boolean list to store the value of checkboxes.
        _checked_list = []
        for check in self.checks:
            _checked_list.append(check.value)

        return _checked_list

    def create_widgets(self):
        """Create widget components in GUI and register callbacks for widgets."""
        self.dpd_shell_no = ipw.Dropdown(
            options=list(range(1, self.no_of_shells + 1)),
            description="Shell No. ",
            value=1,
            layout=ipw.Layout(width="160px"),
        )
        self.dpd_shell_no.observe(self.dpd_shell_no_eventhandler, "value")
        self.btn_prev = ipw.Button(
            icon="chevron-left",
            disabled=True,
            layout=ipw.Layout(width="30px", height="30px"),
        )
        self.btn_prev.on_click(self.on_btn_prev)
        self.btn_next = ipw.Button(
            icon="chevron-right", layout=ipw.Layout(width="30px", height="30px")
        )
        self.btn_next.on_click(self.on_btn_next)

        self.checks = [
            ipw.Checkbox(
                indent=False,
                icon="lock",
                layout=ipw.Layout(
                    width="30px",
                ),
            )
            for element in self.data.elements
        ]
        self.input_items = [
            ipw.BoundedFloatText(min=0, max=1, step=0.01, description=element)
            for element in self.data.elements
        ]
        for i in range(self.no_of_elements):
            self.input_items[i].observe(self.input_item_eventhandler, "value")
            self.input_items[i].index = i
            self.checks[i].observe(self.check_eventhandler, "value")
            self.checks[i].index = i

        self.btn_norm = ipw.Button(
            description="Normalize",
            icon="cog",
            layout=ipw.Layout(width="100px", margin="0 0 0 50px"),
        )
        self.btn_norm.on_click(self.on_btn_norm)
        self.norm_warning = ipw.Valid(
            value=False,
            readout="Unnormalized",
            style={"description_width": "initial"},
            layout=ipw.Layout(visibility="hidden"),
        )

        self.symb_warning = ipw.Valid(
            value=False, layout=ipw.Layout(visibility="hidden")
        )
        self.input_symb = ipw.Text(
            description="Element: ",
            style={"description_width": "initial"},
            placeholder="symbol",
            layout=ipw.Layout(width="125px"),
        )
        self.input_symb.observe(self.input_symb_eventhandler, "value")
        self.btn_add_element = ipw.Button(
            icon="plus-square",
            description="Add",
            disabled=True,
            layout=ipw.Layout(width="60px"),
        )
        self.btn_add_element.on_click(self.on_btn_add_element)

        self.irs_shell_range = ipw.IntRangeSlider(
            min=1,
            max=self.no_of_shells,
            step=1,
            description="Shell No. ",
            disabled=True,
            style={"description_width": "initial"},
            continuous_update=False,
            layout=ipw.Layout(margin="5px 0 0 0"),
        )
        self.irs_shell_range.observe(self.irs_shell_range_eventhandler, "value")

        self.btn_add_shell = ipw.Button(
            icon="plus-square",
            description="Add",
            disabled=True,
            layout=ipw.Layout(
                width="80px",
            ),
        )
        self.btn_add_shell.on_click(self.on_btn_add_shell)
        self.input_v_start = ipw.FloatText(
            min=0,
            description="Add shell(s) with velocity range (km/s): ",
            style={"description_width": "initial"},
            layout=ipw.Layout(
                width="320px",
            ),
        )
        self.input_v_end = ipw.FloatText(
            min=0,
            description="to",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="110px"),
        )
        self.input_v_start.observe(self.input_v_eventhandler, "value")
        self.input_v_end.observe(self.input_v_eventhandler, "value")
        self.overwrite_warning = ipw.HTML(
            value="<font color=darkred><b>Warning:</b></font> Existing shell(s) will be overwritten!",
            layout=ipw.Layout(visibility="hidden", margin="0 0 0 10px"),
        )

        self.btn_output = ipw.Button(
            description="Output CSVY File",
        )
        self.btn_output.on_click(self.on_btn_output)

        self.input_path = ipw.Text(
            description="File path: ", placeholder="Input file name or path"
        )

        self.input_i_time_0 = ipw.FloatText(
            description="Isotope time_0 (day): ",
            style={"description_width": "initial"},
        )

        self.ckb_overwrite = ipw.Checkbox(
            description="overwrite",
            indent=False,
        )

        self.tbs_scale = ipw.ToggleButtons(
            options=["Linear", "Log"],
            description="Scale of yaxes: ",
            style={"description_width": "initial"},
            value="Linear",
        )
        self.tbs_scale.observe(self.tbs_scale_eventhandler, "value")

        self.rbs_single_apply = ipw.RadioButtons(
            options=["Only selected shell"],
            layout=ipw.Layout(margin="10px 0 0 0"),
        )
        self.rbs_single_apply.observe(
            self.rbs_single_apply_eventhandler, "value"
        )
        self.rbs_multi_apply = ipw.RadioButtons(
            options=["A range of shells: "],
            index=None,
            layout=ipw.Layout(width="130px", margin="10px 0 10px 0"),
        )
        self.rbs_multi_apply.observe(self.rbs_multi_apply_eventhandler, "value")

    def update_input_item_value(self, index, value):
        """Update the value of the widget in the list of abundance inputs.

        Keep two decimal places for displayed value and disable
        `input_item_eventhandler` while changing the value.

        Parameters
        ----------
        index : int
            The index of the widget in the list of abundance inputs.
        value : float
            New abundance value.
        """
        self._trigger = False
        # `input_items` is the list of  abundance input widgets.
        self.input_items[index].value = float("{:.2e}".format(value))
        self._trigger = True

    def read_abundance(self):
        """Read abundances data in DataFrame to input items box when
        shell No. changes.
        """
        for i in range(self.no_of_elements):
            value = self.data.abundance.iloc[i, self.shell_no - 1]
            self.update_input_item_value(i, value)

    def bound_locked_sum_to_1(self, index):
        """Ensure the sum of locked abundances is no more than 1. If the
        locked sum is more than 1, calculate the maximum with the sum no
        more than 1 and return it.

        Parameters
        ----------
        index : int
            The index of the widget in the list of abundance inputs.
        """
        locked_mask = np.array(self.checked_list)
        back_value = self.data.abundance.iloc[
            index, self.shell_no - 1
        ]  # abundance value in back end (DataFrame)
        front_value = self.input_items[
            index
        ].value  # abundance value in front end (widget)
        locked_sum = (
            self.data.abundance.loc[locked_mask, self.shell_no - 1].sum()
            - back_value
            + front_value
        )

        if locked_sum > 1:
            new = 1 - (locked_sum - front_value)
            self.data.abundance.iloc[index, self.shell_no - 1] = new
            self.update_input_item_value(index, new)
            self.update_abundance_plot(index)

    def update_abundance_plot(self, index):
        """Update the abundance line of certain element in the plot.

        Parameters
        ----------
        index : int
            The index of the widget in the list of abundance inputs.
        """
        y = self.data.abundance.iloc[index]
        self.fig.data[index + 2].y = np.append(y, y.iloc[-1])

    def update_front_end(self):
        """Update checkbox widgets, input widgets and plot in the front
        end when selected shell No. is changed.
        """
        # Update checkboxes, input boxes and plot.
        for check in self.checks:
            check.value = False
        self.read_abundance()
        self.density_editor.read_density()
        with self.fig.batch_update():
            # Change bar diagonal
            x = list(self.fig.data[0].x)
            width = list(self.fig.data[0].width)
            x_inner = self.data.velocity[self.shell_no - 1].value
            x_outer = self.data.velocity[self.shell_no].value
            x[0] = (x_inner + x_outer) / 2
            self.fig.data[0].x = x
            width[0] = x_outer - x_inner
            self.fig.data[0].width = width

    def update_bar_diagonal(self):
        """Update bar diagonal (and shell no dropdown) when the range
        of shells changed.
        """
        if self.irs_shell_range.disabled:
            x = [self.fig.data[0].x[0]]
            width = [self.fig.data[0].width[0]]
            y = [1]
        else:
            self.shell_no = self.irs_shell_range.value[0]
            self.btn_next.disabled = True
            self.btn_prev.disabled = True

            (start_shell_no, end_shell_no) = self.irs_shell_range.value
            x_inner = self.data.velocity[start_shell_no - 1].value
            x_outer = self.data.velocity[end_shell_no].value
            x = [self.fig.data[0].x[0], (x_outer + x_inner) / 2]
            width = [self.fig.data[0].width[0], x_outer - x_inner]
            y = [1, 1]

        with self.fig.batch_update():
            self.fig.data[0].x = x
            self.fig.data[0].width = width
            self.fig.data[0].y = y

    def update_line_color(self):
        """Update line color in the plot according to colormap."""
        colorscale = transition_colors(self.no_of_elements, self.plot_cmap)
        for i in range(self.no_of_elements):
            self.fig.data[2 + i].line.color = colorscale[i]

    def overwrite_existing_shells(self, v_0, v_1):
        """Judge whether the existing shell(s) will be overwritten when
        inserting a new shell within the entered velocity range.

        Parameters
        ----------
        v_0 : float
            The velocity of inner boundary.
        v_1 : float
            The velocity of outer boundary.

        Returns
        -------
        bool
            True if the existing shell will be overwritten, False otherwise.
        """
        v_vals = self.data.velocity.value
        position_0 = np.searchsorted(v_vals, v_0)
        position_1 = np.searchsorted(v_vals, v_1)

        index_0 = (
            position_0 - 1
            if np.isclose(v_vals[position_0 - 1], v_0)
            else position_0
        )
        index_1 = (
            position_1 - 1
            if np.isclose(v_vals[position_1 - 1], v_1)
            else position_1
        )

        if (index_1 - index_0 > 1) or (
            index_1 - index_0 == 1 and not np.isclose(v_vals[index_0], v_0)
        ):
            return True
        else:
            return False

    def on_btn_add_shell(self, obj):
        """Add new shell with given boundary velocities. Triggered if
        the button is clicked.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        v_start = self.input_v_start.value
        v_end = self.input_v_end.value
        v_vals = self.data.velocity.value

        position_0 = np.searchsorted(v_vals, v_start)
        position_1 = np.searchsorted(v_vals, v_end)
        start_index = (
            int(position_0 - 1)
            if np.isclose(v_vals[position_0 - 1], v_start)
            else int(position_0)
        )
        end_index = (
            int(position_1 - 1)
            if np.isclose(v_vals[position_1 - 1], v_end)
            else int(position_1)
        )

        if (
            end_index < self.no_of_shells
            and np.isclose(v_vals[start_index], v_start)
            and np.isclose(v_vals[end_index], v_end)
            and end_index - start_index == 1
        ):
            return

        if start_index < self.no_of_shells and np.isclose(
            v_vals[start_index], v_start
        ):
            new_shell_abundances = self.data.abundance[start_index]
        else:
            index = min(start_index - 1, self.no_of_shells - 1)
            new_shell_abundances = self.data.abundance[max(index, 0)]

        # Delete the overwritten shell (abundances and velocities).
        if end_index < len(v_vals) and np.isclose(v_vals[end_index], v_end):
            # New shell will overwrite the original shell that ends at v_end.
            v_vals = np.delete(v_vals, end_index)
            self.data.abundance.drop(max(0, end_index - 1), 1, inplace=True)

        # Insert new velocities calculate new densities according
        # to new velocities through interpolation.
        v_vals = np.insert(v_vals, [start_index, end_index], [v_start, v_end])
        v_vals = np.delete(v_vals, slice(start_index + 1, end_index + 1))

        self.data.density = (
            np.interp(
                v_vals,
                self.data.velocity[1:].value,
                self.data.density[1:].value,
            )
            * self.data.density.unit
        )
        self.data.velocity = v_vals * self.data.velocity.unit

        # Change abundances after adding new shell.
        if start_index != end_index:
            self.data.abundance.insert(start_index, "", new_shell_abundances)
            self.data.abundance.drop(
                self.data.abundance.iloc[:, start_index : end_index - 1],
                1,
                inplace=True,
            )
        else:
            if start_index == 0:
                self.data.abundance.insert(0, "new", new_shell_abundances)
                self.data.abundance.insert(
                    0, "gap", new_shell_abundances
                )  # Add a shell to fill the gap.
            else:
                self.data.abundance.insert(
                    start_index - 1, "new", new_shell_abundances
                )
                self.data.abundance.insert(
                    start_index - 1, "gap", new_shell_abundances
                )  # Add a shell to fill the gap with original abundances

        self.data.abundance.columns = range(self.no_of_shells)

        # Update data and x axis in plot.
        with self.fig.batch_update():
            self.fig.layout.xaxis.autorange = True
            self.fig.data[1].x = self.data.velocity
            self.fig.data[1].y = np.append(
                self.data.density[1:], self.data.density[-1]
            )
            for i in range(self.no_of_elements):
                self.fig.data[i + 2].x = self.data.velocity
                self.update_abundance_plot(i)

        self.dpd_shell_no.options = list(range(1, self.no_of_shells + 1))
        self.update_bar_diagonal()
        self.update_front_end()
        self.irs_shell_range.max = self.no_of_shells
        self.overwrite_warning.layout.visibility = "hidden"

    def tbs_scale_eventhandler(self, obj):
        """Switch the scale type of y axis between linear mode and log
        mode. Triggered if the toggle button is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        scale_mode = obj.new

        if scale_mode == "Linear":
            self.fig.update_layout(
                yaxis=dict(
                    type="linear",
                    range=[0, 1],
                ),
                yaxis2=dict(type="linear"),
            )
        else:
            self.fig.update_layout(
                yaxis=dict(
                    type="log",
                    range=[-8, 0],
                ),
                yaxis2=dict(type="log"),
            )

    def input_item_eventhandler(self, obj):
        """Update the data and the widget when it gets new abundance
        input. Triggered if the abundance input is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        if self._trigger:
            item_index = obj.owner.index
            is_locked = self.checks[item_index]

            if is_locked:
                self.bound_locked_sum_to_1(item_index)

            self.data.abundance.iloc[
                item_index, self.shell_no - 1
            ] = obj.owner.value

            if self.rbs_multi_apply.index is None:
                self.update_abundance_plot(item_index)
            else:
                self.apply_to_multiple_shells(item_index)

        if np.isclose(self.data.abundance.iloc[:, self.shell_no - 1].sum(), 1):
            self.norm_warning.layout.visibility = "hidden"
        else:
            self.norm_warning.layout.visibility = "visible"

    def check_eventhandler(self, obj):
        """Triggered if the checkbox is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        item_index = obj.owner.index

        if obj.new == True:
            self.bound_locked_sum_to_1(item_index)

    def dpd_shell_no_eventhandler(self, obj):
        """Make the data in widgets correspond with the selected shell.
        Triggered if the dropdown value is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        # Disable "previous" and "next" buttons when shell no comes to boundaries.
        if obj.new == 1:
            self.btn_prev.disabled = True
        else:
            self.btn_prev.disabled = False

        if obj.new == self.no_of_shells:
            self.btn_next.disabled = True
        else:
            self.btn_next.disabled = False

        self.update_front_end()

    def on_btn_prev(self, obj):
        """Move to previous shell.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        self.shell_no -= 1

    def on_btn_next(self, obj):
        """Move to next shell.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        self.shell_no += 1

    def on_btn_norm(self, obj):
        """Normalize unlocked abundances to 1. Triggered if the
        normalize button is clicked.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        locked_mask = np.array(self.checked_list)
        locked_sum = self.data.abundance.loc[
            locked_mask, self.shell_no - 1
        ].sum()
        unlocked_arr = self.data.abundance.loc[~locked_mask, self.shell_no - 1]

        # if abundances are all zero
        if unlocked_arr.sum() == 0:
            return

        self.data.abundance.loc[~locked_mask, self.shell_no - 1] = (
            (1 - locked_sum) * unlocked_arr / unlocked_arr.sum()
        )

        self.read_abundance()

        for i in range(self.no_of_elements):
            if self.rbs_multi_apply.index is None:
                # Only normalize selected shell
                self.update_abundance_plot(i)
            else:
                self.apply_to_multiple_shells(i)

    @debounce(0.5)
    def input_symb_eventhandler(self, obj):
        """Judge whether the input symbol is valid. Triggered after 0.5s
        when the symbol input is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        element_symbol_string = obj.new.capitalize()

        if element_symbol_string == "":
            self.symb_warning.layout.visibility = "hidden"
            self.btn_add_element.disabled = True
            return

        if element_symbol_string in self.data.elements:
            self.symb_warning.layout.visibility = "visible"
            self.btn_add_element.disabled = True
            self.symb_warning.readout = "Already exists!"
            return

        try:
            if is_valid_nuclide_or_elem(element_symbol_string):
                self.symb_warning.layout.visibility = "hidden"
                self.btn_add_element.disabled = False
                return

        except RuntimeError:
            pass

        self.symb_warning.layout.visibility = "visible"
        self.btn_add_element.disabled = True
        self.symb_warning.readout = "invalid"

    def on_btn_add_element(self, obj):
        """Add new element and update the display in the front end.
        Triggered if the add button is clicked.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        element_symbol_string = self.input_symb.value.capitalize()

        if element_symbol_string in Z_DICT.values():
            z = elem_to_Z(element_symbol_string)
            self.data.abundance.loc[(z, ""), :] = 0
        else:
            nuc = Nuclide(element_symbol_string)
            mass_no = nuc.A
            z = nuc.Z
            self.data.abundance.loc[(z, mass_no), :] = 0

        self.data.abundance.sort_index(inplace=True)

        # Add new BoundedFloatText control and Checkbox control.
        item = ipw.BoundedFloatText(min=0, max=1, step=0.01)
        check = ipw.Checkbox(indent=False, layout=ipw.Layout(width="30px"))
        item.index = self.no_of_elements - 1
        check.index = self.no_of_elements - 1
        item.observe(self.input_item_eventhandler, "value")
        check.observe(self.check_eventhandler, "value")
        self.input_items.append(item)
        self.checks.append(check)

        # Keep the order of description same with atomic number
        self.data.elements = self.data.get_symbols()
        for i in range(self.no_of_elements):
            self.input_items[i].description = self.data.elements[i]

        self.box_editor.children = [
            ipw.VBox(self.input_items),
            ipw.VBox(self.checks),
        ]

        with self.fig.batch_update():
            # Add new trace to plot.
            self.fig.add_scatter(
                x=self.data.velocity,  # convert to km/s
                y=[0] * (self.no_of_shells + 1),
                mode="lines+markers",
                name=element_symbol_string,
                line=dict(shape="hv"),
            )
            # Sort the legend in atomic order.
            fig_data_lst = list(self.fig.data)
            fig_data_lst.insert(
                np.argwhere(self.data.elements == element_symbol_string)[0][0]
                + 2,
                self.fig.data[-1],
            )
            self.fig.data = fig_data_lst[:-1]

            self.update_line_color()

        self.read_abundance()

        # Clear symbol input box.
        self.input_symb.value = ""

    def apply_to_multiple_shells(self, item_index):
        """Apply the changed abundances to specified range of shell(s).

        Parameters
        ----------
        item_index : int
            The index of the widget in the list of abundance inputs.
        """
        start_index = self.irs_shell_range.value[0] - 1
        end_index = self.irs_shell_range.value[1]
        applied_shell_index = self.shell_no - 1

        self.data.abundance.iloc[
            item_index, start_index:end_index
        ] = self.data.abundance.iloc[item_index, applied_shell_index]

        self.update_abundance_plot(item_index)

    def input_v_eventhandler(self, obj):
        """Judge whether the input velocity range is valid. Triggered if the
        velocity input is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        v_start = self.input_v_start.value
        v_end = self.input_v_end.value

        if v_start >= v_end or v_start < 0 or v_end < 0:
            self.btn_add_shell.disabled = True
            return
        else:
            self.btn_add_shell.disabled = False

        # Whether overwrite existing shells
        if self.overwrite_existing_shells(v_start, v_end):
            self.overwrite_warning.layout.visibility = "visible"
        else:
            self.overwrite_warning.layout.visibility = "hidden"

    def on_btn_output(self, obj):
        """Triggered if the output button is clicked.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        path = self.input_path.value
        overwrite = self.ckb_overwrite.value

        self.to_csvy(path, overwrite)

    def rbs_single_apply_eventhandler(self, obj):
        """Switch to single shell editing mode. Triggered if the
        first radio button is selected.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        self.abundance_note.layout.visibility = "hidden"

        self.rbs_multi_apply.unobserve(
            self.rbs_multi_apply_eventhandler, "value"
        )
        self.rbs_multi_apply.index = None
        self.rbs_multi_apply.observe(self.rbs_multi_apply_eventhandler, "value")

        self.dpd_shell_no.disabled = False
        self.btn_next.disabled = False
        self.btn_prev.disabled = False
        self.irs_shell_range.disabled = True
        self.update_bar_diagonal()

    def rbs_multi_apply_eventhandler(self, obj):
        """Switch to multi-shells editing mode. Triggered if the
        second radio button is selected.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        self.abundance_note.layout.visibility = "visible"

        self.shell_no = self.irs_shell_range.value[0]

        self.rbs_single_apply.unobserve(
            self.rbs_single_apply_eventhandler, "value"
        )
        self.rbs_single_apply.index = None
        self.rbs_single_apply.observe(
            self.rbs_single_apply_eventhandler, "value"
        )
        # self.irs_shell_range.value = [self.shell_no, self.shell_no]
        self.dpd_shell_no.disabled = True
        self.btn_next.disabled = True
        self.btn_prev.disabled = True
        self.irs_shell_range.disabled = False
        self.update_bar_diagonal()

    def irs_shell_range_eventhandler(self, obj):
        """Select the velocity range of new shell and highlight the range
        in the plot. Triggered if the shell range slider is changed.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        self.update_bar_diagonal()

    def generate_abundance_density_plot(self):
        """Generate abundance and density plot in different shells."""
        title = "Abundance/Density vs Velocity"
        abundance = self.data.abundance
        velocity = self.data.velocity
        density = self.data.density

        # Bar Diagonal
        self.fig.add_trace(
            go.Bar(
                x=[(velocity[0].value + velocity[1].value) / 2],
                y=[1],
                width=[velocity[1].value - velocity[0].value],
                name="Selected shell",
                marker=dict(
                    color="rgb(253,205,172)",
                ),
                hoverinfo="none",
            )
        )

        self.fig.add_trace(
            go.Scatter(
                x=velocity,
                y=np.append(density[1:], density[-1]),
                mode="lines+markers",
                name="<b>Density</b>",
                yaxis="y2",
                line=dict(color="black", shape="hv"),
                marker_symbol="square",
            ),
        )

        for i in range(self.no_of_elements):
            self.fig.add_trace(
                go.Scatter(
                    x=velocity,
                    y=np.append(abundance.iloc[i], abundance.iloc[i, -1]),
                    mode="lines+markers",
                    line=dict(shape="hv"),
                    name=self.data.elements[i],
                ),
            )

        self.fig.update_layout(
            xaxis=dict(
                title="Velocity (km/s)",
                tickformat="f",
            ),
            yaxis=dict(
                title="Fractional Abundance",
                exponentformat="e",
                range=[0, 1],
            ),
            yaxis2=dict(
                title="Density (g/cm^3)",
                exponentformat="e",
                overlaying="y",
                side="right",
            ),
            height=500,
            title=title,
            hovermode="closest",
            legend=dict(
                x=1.15,
            ),
        )

    def display(self, cmap="jet"):
        """Display the GUI.

        Parameters
        ----------
        cmap : str, default: "jet", optional
            String defines the colormap used in abundance density plot.

        Returns
        -------
        ipywidgets.widgets.widget_box.VBox
            A box that contains all the widgets in the GUI.
        """
        if not is_notebook():
            print("Please use a notebook to display the widget")
        else:
            # --------------Combine widget components--------------
            self.box_editor = ipw.HBox(
                [
                    ipw.VBox(self.input_items),
                    ipw.VBox(
                        self.checks, layout=ipw.Layout(margin="0 0 0 10px")
                    ),
                ]
            )

            box_add_shell = ipw.HBox(
                [
                    self.input_v_start,
                    self.input_v_end,
                    self.btn_add_shell,
                    self.overwrite_warning,
                ],
                layout=ipw.Layout(margin="0 0 0 50px"),
            )

            box_head = ipw.HBox(
                [self.dpd_shell_no, self.btn_prev, self.btn_next, box_add_shell]
            )

            box_add_element = ipw.HBox(
                [self.input_symb, self.btn_add_element, self.symb_warning],
                layout=ipw.Layout(margin="0 0 0 80px"),
            )

            help_note = ipw.HTML(
                value="<p style='text-indent: 40px'>* Select a checkbox "
                "to lock the abundance of corresponding element. </p>"
                "<p style='text-indent: 40px'> On clicking the 'Normalize' "
                "button, the locked abundance(s) will <b>not be normalized</b>."
                " </p>",
                indent=True,
            )

            self.abundance_note = ipw.HTML(
                description="(The following abundances are for the innermost "
                "shell in selected range.)",
                layout=ipw.Layout(visibility="hidden"),
                style={"description_width": "initial"},
            )

            box_norm = ipw.HBox([self.btn_norm, self.norm_warning])

            box_apply = ipw.VBox(
                [
                    ipw.Label(value="Apply abundance(s) to:"),
                    self.rbs_single_apply,
                    ipw.HBox(
                        [
                            self.rbs_multi_apply,
                            self.irs_shell_range,
                            self.abundance_note,
                        ]
                    ),
                ],
                layout=ipw.Layout(margin="0 0 15px 50px"),
            )

            box_features = ipw.VBox([box_norm, help_note])
            box_abundance = ipw.VBox(
                [
                    box_apply,
                    ipw.HBox([self.box_editor, box_features]),
                    box_add_element,
                ]
            )
            box_density = self.density_editor.display()

            main_tab = ipw.Tab([box_abundance, box_density])
            main_tab.set_title(0, "Edit Abundance")
            main_tab.set_title(1, "Edit Density")

            hint = ipw.HTML(
                value="<b><font size='3'>Save model as file: </font></b>"
            )
            box_output = ipw.VBox(
                [
                    hint,
                    self.input_i_time_0,
                    ipw.HBox(
                        [self.input_path, self.btn_output, self.ckb_overwrite]
                    ),
                ]
            )

            # Initialize the widget and plot colormap
            self.plot_cmap = cmap
            self.update_line_color()
            self.read_abundance()
            self.density_editor.read_density()

            return ipw.VBox(
                [
                    self.tbs_scale,
                    self.fig,
                    box_head,
                    main_tab,
                    box_output,
                    self.error_view,
                ]
            )

    @error_view.capture(clear_output=True)
    def to_csvy(self, path, overwrite):
        """Output CSVY file on the specified path.

        Parameters
        ----------
        path : str
            Output path.
        overwrite : bool
            True if overwriting, False otherwise.
        """
        posix_path = Path(path)
        posix_path = posix_path.with_suffix(".csvy")

        if posix_path.exists() and not overwrite:
            raise FileExistsError(
                "The file already exists. Click the 'overwrite' checkbox to overwrite it."
            )
        else:
            self.write_yaml_portion(posix_path)
            self.write_csv_portion(posix_path)

    @error_view.capture(clear_output=True)
    def write_yaml_portion(self, path):
        """Write the YAML portion of the output file.

        Parameters
        ----------
        path : pathlib.PosixPath
        """
        name = path.name
        d_time_0 = self.data.density_t_0
        i_time_0 = self.input_i_time_0.value * u.day
        custom_yaml = CustomYAML(
            name,
            d_time_0,
            i_time_0,
            self.data.velocity[0],
            self.data.velocity[-1],
        )
        custom_yaml.create_fields_dict(self.data.elements)

        with path.open("w") as f:
            yaml_output = yaml.dump(custom_yaml, sort_keys=False)

            # Add YAML delimiter
            yaml_output_snippets = yaml_output.split("\n")
            yaml_output_snippets[0] = yaml_output_snippets[-1] = YAML_DELIMITER
            yaml_output = "\n".join(yaml_output_snippets) + "\n"

            f.write(yaml_output)

        print("Saved Successfully!")

    def write_csv_portion(self, path):
        """Write the CSV portion of the output file.

        Parameters
        ----------
        path : pathlib.PosixPath
        """
        try:
            data = self.data.abundance.T
            data.columns = self.data.elements
            first_row = [0] * self.no_of_elements
            data.loc[-1] = first_row
            data.index += 1  # shifting index
            data.sort_index(inplace=True)

            formatted_v = pd.Series(self.data.velocity.value).apply(
                lambda x: "%.3e" % x
            )
            # Make sure velocity is within the boundary.
            formatted_v[0] = self.data.velocity.value[0]
            formatted_v[-1] = self.data.velocity.value[-1]

            density = self.data.density
            data.insert(0, "velocity", formatted_v)
            data.insert(1, "density", density)

            data.to_csv(path, mode="a", index=False)
        except pd.errors.EmptyDataError:
            data = None

    @classmethod
    def from_csvy(cls, fpath):
        """Create a new CustomAbundanceWidget instance with data from CSVY file.

        Parameters
        ----------
        fpath : str
            the path of CSVY file.

        Returns
        -------
        CustomAbundanceWidget
        """
        widget_data = CustomAbundanceWidgetData.from_csvy(fpath)

        return cls(widget_data)

    @classmethod
    def from_yml(cls, fpath):
        """Create a new CustomAbundanceWidget instance with data from YAML file.

        Parameters
        ----------
        fpath : str
            The path of YAML file.

        Returns
        -------
        CustomAbundanceWidget
        """
        widget_data = CustomAbundanceWidgetData.from_yml(fpath)

        return cls(widget_data)

    @classmethod
    def from_hdf(cls, fpath):
        """Create a new CustomAbundanceWidget instance with data from HDF file.

        Parameters
        ----------
        fpath : str
            the path of HDF file.

        Returns
        -------
        CustomAbundanceWidget
        """
        widget_data = CustomAbundanceWidgetData.from_hdf(fpath)

        return cls(widget_data)

    @classmethod
    def from_simulation(cls, sim):
        """Create a new CustomAbundanceWidget instance from a Simulation object.

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        CustomAbundanceWidget
        """
        widget_data = CustomAbundanceWidgetData.from_simulation(sim)

        return cls(widget_data)


class DensityEditor:
    """Widget to edit density profile of the model.

    It provides an interface to allow the user directly change
    the density, or calculate the density with given type and
    parameters.

    Attributes
    ----------
    shell_no : int
        The selected shell number.
    _trigger : bool
        If False, disable the callback when density input is changed.
    """

    def __init__(self, widget_data, fig, shell_no_widget):
        """Initialize DensityEditor with data and widget components.

        Parameters
        ----------
        widget_data : CustomAbundanceWidgetData
            Data in the custom abundance widget.
        fig : plotly.graph_objs._figurewidget.FigureWidget
            The figure object of density plot.
        shell_no_widget : ipywidgets.widgets.widget_selection.Dropdown
            A widget to record the selected shell number.
        """
        self.data = widget_data
        self.fig = fig
        self.shell_no_widget = shell_no_widget

        self.create_widgets()
        self._trigger = True

    @property
    def shell_no(self):
        return self.shell_no_widget.value

    def create_widgets(self):
        """Create widget components in density editor GUI and register
        callbacks for widgets.
        """
        self.input_d = ipw.FloatText(
            description="Density",
            layout=ipw.Layout(
                width="230px",
            ),
        )
        self.input_d.observe(self.input_d_eventhandler, "value")

        self.input_d_time_0 = ipw.FloatText(
            value=self.data.density_t_0.value,
            description="Density time_0 (day): ",
            style={"description_width": "initial"},
            layout=ipw.Layout(margin="0 0 20px 0"),
        )
        self.input_d_time_0.observe(self.input_d_time_0_eventhandler, "value")

        self.input_d_time_0 = ipw.FloatText(
            value=self.data.density_t_0.value,
            description="Density time_0 (day): ",
            style={"description_width": "initial"},
            layout=ipw.Layout(margin="0 0 20px 0"),
        )
        self.input_d_time_0.observe(self.input_d_time_0_eventhandler, "value")

        self.dpd_dtype = ipw.Dropdown(
            options=["-", "uniform", "exponential", "power_law"],
            description="Density type: ",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="300px"),
        )
        self.dpd_dtype.observe(self.dpd_dtype_eventhandler, "value")

        self.input_rho_0 = ipw.FloatText(
            description="rho_0", layout=ipw.Layout(width="300px")
        )

        self.input_exp = ipw.FloatText(
            description="exponent", layout=ipw.Layout(width="300px")
        )

        self.input_v_0 = ipw.FloatText(
            description="v_0", layout=ipw.Layout(width="300px")
        )

        self.input_value = ipw.FloatText(
            description="value", layout=ipw.Layout(width="300px")
        )

        self.btn_calculate = ipw.Button(
            icon="calculator",
            description="Calculate Density",
            layout=ipw.Layout(width="300px"),
        )
        self.btn_calculate.on_click(self.on_btn_calculate)

        self.uniform_box = ipw.HBox(
            [self.input_value, ipw.Label(value="g cm^3")]
        )

        # Formula to compute density profile
        form_exp = ipw.HTMLMath(
            value=r"$$\rho = \rho_0 \times \exp \left( -\frac{v}{v_0} \right)$$",
            layout=ipw.Layout(margin="0 0 0 100px"),
        )
        form_pow = ipw.HTMLMath(
            value=r"$$\rho = \rho_0 \times \left( \frac{v}{v_0} \right)^n$$",
            layout=ipw.Layout(margin="0 0 0 100px"),
        )

        self.exp_box = ipw.VBox(
            [
                form_exp,
                ipw.HBox([self.input_rho_0, ipw.Label(value="g cm^3")]),
                ipw.HBox([self.input_v_0, ipw.Label(value="km/s")]),
            ]
        )
        self.pow_box = ipw.VBox(
            [
                form_pow,
                ipw.HBox([self.input_rho_0, ipw.Label(value="g cm^3")]),
                ipw.HBox([self.input_v_0, ipw.Label(value="km/s")]),
                self.input_exp,
            ]
        )

    def read_density(self):
        """Read density data in DataFrame to density input box when
        shell No. changes.
        """
        dvalue = self.data.density[self.shell_no].value
        self._trigger = False
        self.input_d.value = float("{:.3e}".format(dvalue))
        self._trigger = True

    def update_density_plot(self):
        """Update the density line in the plot."""
        y = np.append(self.data.density[1:], self.data.density[-1])
        self.fig.data[1].y = y

    def input_d_eventhandler(self, obj):
        """Update the data and the widgets when the widget gets new
        density input.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        if self._trigger:
            new_value = obj.new
            self.data.density[self.shell_no] = (
                new_value * self.data.density.unit
            )

            self.update_density_plot()

    def input_d_time_0_eventhandler(self, obj):
        """Update density time 0 data when the widget gets new input.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        new_value = obj.new
        self.data.density_t_0 = new_value * self.data.density_t_0.unit

    def input_d_time_0_eventhandler(self, obj):
        """Update density time 0 data when the widget gets new input.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.
        """
        new_value = obj.new
        self.data.density_t_0 = new_value * self.data.density_t_0.unit

    dtype_out = ipw.Output()

    @dtype_out.capture(clear_output=True)
    def dpd_dtype_eventhandler(self, obj):
        """Display the input boxes of the specified density type.
        Triggered if the density type dropdown is changed.

        Parameters
        ----------
        obj : traitlets.utils.bunch.Bunch
            A dictionary holding the information about the change.

        Returns
        -------
        ipywidgets.widgets.widget_box.VBox
            A box widget that contains the input boxes of certain density
            type parameters.
        """
        if obj.new == "uniform":
            display(self.uniform_box)
        elif obj.new == "exponential":
            display(self.exp_box)
        elif obj.new == "power_law":
            display(self.pow_box)

    def on_btn_calculate(self, obj):
        """Calculate density according to density parameters input.
        Triggered if the density calculate button is clicked.

        Parameters
        ----------
        obj : ipywidgets.widgets.widget_button.Button
            The clicked button instance.
        """
        dtype = self.dpd_dtype.value
        density = self.data.density
        velocity = self.data.velocity

        if dtype == "-":
            return

        if dtype == "uniform":
            if self.input_value.value == 0:
                return

            self.data.density = (
                self.input_value.value * density.unit * np.ones(len(density))
            )
        else:
            if self.input_v_0.value == 0 or self.input_rho_0.value == 0:
                return

            adjusted_velocity = velocity.insert(0, 0)
            v_middle = (
                adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
            )
            v_0 = self.input_v_0.value * velocity.unit
            rho_0 = self.input_rho_0.value * density.unit

            if dtype == "exponential":
                self.data.density = calculate_exponential_density(
                    v_middle, v_0, rho_0
                )

            elif dtype == "power_law":
                exponent = self.input_exp.value
                self.data.density = calculate_power_law_density(
                    v_middle, v_0, rho_0, exponent
                )

        self.read_density()
        self.update_density_plot()

    def display(self):
        """Display the GUI.

        Returns
        -------
        ipywidgets.widgets.widget_box.VBox
            A box that contains all the widgets in the GUI.
        """
        hint1 = ipw.HTML(
            value="<font size='3'>1) Edit density of the selected shell:</font>"
        )
        hint2 = ipw.HTML(
            value="<font size='3'>2) Edit densities for all shells:</font>"
        )
        d_box = ipw.HBox(
            [self.input_d, ipw.Label(value="g/cm^3")],
            layout=ipw.Layout(margin="0 0 20px 0"),
        )
        return ipw.VBox(
            [
                self.input_d_time_0,
                hint1,
                d_box,
                hint2,
                self.dpd_dtype,
                self.dtype_out,
                self.btn_calculate,
            ]
        )
