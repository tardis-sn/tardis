"""Class to create and display Custom Abundance Widget."""
import os
import math
import numpy as np
import pandas as pd
import ipywidgets as ipw
import plotly.graph_objects as go
from astropy import units as u
from pyne import nucname

from tardis.util.base import quantity_linspace
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.model.density import calculate_power_law_density, calculate_exponential_density
from tardis.io.config_validator import validate_dict
from tardis.io.parsers.csvy import load_csvy
from tardis.io.model_reader import (
    read_uniform_abundances,
    parse_csv_abundances,
)
from tardis.util.base import (
    atomic_number2element_symbol,
    quantity_linspace
)

import asyncio

class Timer:
    """
    Cited from https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html
    """
    def __init__(self, timeout, callback):
        self._timeout = timeout
        self._callback = callback

    async def _job(self):
        await asyncio.sleep(self._timeout)
        self._callback()

    def start(self):
        self._task = asyncio.ensure_future(self._job())

    def cancel(self):
        self._task.cancel()

def debounce(wait):
    """ 
    Decorator that will postpone a function's
    execution until after `wait` seconds
    have elapsed since the last time it was invoked. 
    Cited from https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html
    """
    def decorator(fn):
        timer = None
        def debounced(*args, **kwargs):
            nonlocal timer
            def call_it():
                fn(*args, **kwargs)
            if timer is not None:
                timer.cancel()
            timer = Timer(wait, call_it)
            timer.start()
        return debounced
    return decorator

class CustomAbundanceWidget:
    out_plot = ipw.Output()
    def __init__(self, density, abundance, velocity):
        self.density = density.to("g cm^-3")
        self.abundance = abundance
        self.velocity = velocity.to("km/s")
        self.elements = self.get_symbols()

        self._trigger = True # Whether to trigger 'input_item_eventhandler' when 'input_item' is changed 

        self.create_widgets()

    @property
    def shell_no(self):
        return self.dd_shell_no.value

    @shell_no.setter
    def shell_no(self, value):
        self.dd_shell_no.value = value

    @property
    def no_of_shells(self):
        return self.abundance.shape[1]

    @property
    def no_of_elements(self):
        return self.abundance.shape[0]

    @property
    def checked_list(self): # A boolean list to store the value of checkboxes.
        _checked_list = []
        for check in self.checks:
            _checked_list.append(check.value)

        return _checked_list

    def get_symbols(self):
        str_symbols = np.array(self.abundance.index.get_level_values(0).map(atomic_number2element_symbol))
        str_mass = np.array(self.abundance.index.get_level_values(1), dtype="str")
        return(np.add(str_symbols, str_mass))
        
    def create_widgets(self):
        self.dd_shell_no = ipw.Dropdown(options=list(range(1, self.no_of_shells+1)), 
                                description="Shell No. ", 
                                value=1,
                                layout=ipw.Layout(width="160px")
                            )
        self.dd_shell_no.observe(self.dd_shell_no_eventhandler, "name")

        self.btn_prev = ipw.Button(icon="chevron-left", 
                                disabled=True, 
                                layout=ipw.Layout(width="30px", height="30px")
                                )
        self.btn_next = ipw.Button(icon="chevron-right", 
                            layout=ipw.Layout(width="30px", height="30px")
                            )

        self.checks = [ipw.Checkbox(indent=False, 
                            layout=ipw.Layout(width="30px")
                            ) for element in self.elements]

        self.input_items = [ipw.BoundedFloatText(min=0, 
                                    max=1, 
                                    step=0.01, 
                                    description=element) 
                for element in self.elements]

        for i in self.no_of_elements:
            self.input_items[i].observe(self.input_item_eventhandler, "value")
            self.input_items[i].index = i
            self.checks[i].observe(self.check_handler, "value")
            self.checks[i].index = i

        self.btn_norm = ipw.Button(description="Normalize",
                                icon="cog", 
                                layout = ipw.Layout(width="100px",
                                                    margin="0 0 0 50px")
                                )
        self.btn_norm.on_click(self.on_btn_norm)

        self.tb_scale = ipw.ToggleButtons(
            options=['Linear', 'Log'],
            description='Scale of yaxes: ',
            style={'description_width': 'initial'},
            value='Linear'
        )
        self.tb_scale.observe(self.scale_handler, "value")

        self.norm_warning = ipw.Valid(
            value=False,
            readout='Unnormalized',
            style={'description_width': 'initial'},
            layout=ipw.Layout(visibility='hidden')
        )

        self.symb_warning = ipw.Valid(
            value=False,
            layout=ipw.Layout(visibility='hidden')
        )

        self.input_symb = ipw.Text(
            description='Element: ',
            style={'description_width': 'initial'},
            placeholder='symbol',
            layout=ipw.Layout(width='120px'),
        )
        self.input_symb.observe(self.input_symb_eventhandler, "value")

        self.btn_add_element = ipw.Button(
            icon='plus-square', 
            description='Add',
            disabled=True,
            layout=ipw.Layout(
                width='60px'
            )
        )
        self.btn_add_element.on_click(self.on_btn_add_element)





    def update_input_item_value(self, index, value):
        self._trigger = False
        self.input_items[index].value = float("{:.2e}".format(value))
        self._trigger = True        

    # Read abundances to input items box when shell no changes
    def read_abundance(self):
        for i in range(self.no_of_elements):
            value = self.abundance.iloc[i,self.shell_no-1]
            self.update_input_item_value(i, value)

    # Ensure the sum of locked elements less than 1
    def bound_locked_sum_to_1(self, index):
        locked_mask = np.array(self.checked_list)
        back_value = self.abundance[index, self.shell_no-1] # abundance value in back end (DataFrame)
        front_value = self.input_items[index].value # abundance value in front end (widget)
        locked_sum = self.abundance.loc[locked_mask, self.shell_no-1].sum() - back_value + front_value

        if locked_sum > 1:
            new = 1 - (locked_sum - back_value)
            self.abundance[index, self.shell_no-1] = new
            self.update_input_item_value(self, index, new)

    def scale_eventhandler(self, obj):
        scale_mode = obj.new

        if scale_mode == 'Linear':
            self.fig.update_layout(
                yaxis=dict(
                        type="linear",
                        range=[0, 1],
                        ),
                yaxis2=dict(type="linear")
            )
        else:
            self.fig.update_layout(
                yaxis=dict(
                        type="log",
                        range=[-8, 0],
                        ),
                yaxis2=dict(type="log")
            )

    def input_item_eventhandler(self, obj):
        if self._trigger:
            item_index = obj.owner.index
            is_locked = self.checks[item_index]
            
            if is_locked:
                self.bound_locked_sum_to_1(item_index)

            self.abundance.iloc[item_index, self.shell_no-1] = obj.new
            
            if math.isclose(self.abundance.iloc[:, self.shell_no-1].sum(), 1):
                self.norm_warning.layout.visibility = 'hidden'
            else:
                self.norm_warning.layout.visibility = 'visible'
            
            # Update plot
            self.fig.data[item_index+1].y = self.abundance.iloc[item_index]

    def check_eventhandler(self, obj):
        item_index = obj.owner.index

        if obj.new == True:
            self.bound_locked_sum_to_1(item_index)

    def dd_shell_no_eventhandler(self, obj):  
        # Disable 'previous' and 'next' buttons when shell no comes to boundaries.
        if obj.new == 1:
            self.btn_prev.disabled = True
        else:
            self.btn_prev.disabled = False
        
        if obj.new == self.no_of_shells:
            self.btn_next.disabled = True
        else:
            self.btn_next.disabled = False

        # Update checkboxes, input boxes and plot.
        for check in self.checks:
            check.value = False
        self.read_abundance()
        self.read_density()
        with self.fig.batch_update():
            # Change line diagonal
            self.fig.layout.shapes[0].x0 = self.velocity[self.shell_no].value
            self.fig.layout.shapes[0].x1 = self.velocity[self.shell_no].value  

    def on_btn_prev(self, obj):
        self.shell_no -= 1
    
    def on_btn_next(self, obj):
        self.shell_no += 1

    def on_btn_norm(self, obj):
        locked_mask = np.array(self.checked_list)
        locked_sum = self.abundance.loc[locked_mask, self.shell_no-1].sum()
        unlocked_arr = self.abundance.loc[~locked_mask, self.shell_no-1]
        
        # if abundances are all zero
        if unlocked_arr.sum() == 0:
            return
        
        self.abundance.loc[~locked_mask, self.shell_no-1] = (1 - locked_sum) * unlocked_arr / unlocked_arr.sum()
        
        self.read_abundance()

    @debounce(0.5)
    def input_symb_eventhandler(self, obj):
        element_symbol_string = obj.new.capitalize()        
        
        if element_symbol_string == "":
            self.symb_warning.layout.visibility = 'hidden' 
            self.btn_add_element.disabled = True
            return
        
        try:
            if nucname.iselement(element_symbol_string) or nucname.isnuclide(element_symbol_string):
                self.symb_warning.layout.visibility = 'hidden'
                self.btn_add_element.disabled = False
                return

        except RuntimeError:
            pass        
        
        self.symb_warning.layout.visibility = 'visible'
        self.btn_add_element.disabled = True

        if element_symbol_string in self.elements:        
            self.symb_warning.readout = 'Already exists!'
        else:
            self.symb_warning.readout = 'invalid'
    
    def on_btn_add_element(self, obj):
        element_symbol_string = self.input_symb.value.capitalize()
        
        if element_symbol_string in nucname.name_zz:
            z = nucname.name_zz[element_symbol_string]
            self.abundance.loc[(z, ''), :] = 0
        else:
            mass_no = nucname.anum(element_symbol_string)
            z = nucname.znum(element_symbol_string)
            self.abundance.loc[(z, mass_no), :] = 0
        
        self.abundance.sort_index(inplace=True)
        
        # Add new BoundedFloatText control and Checkbox control.
        _item = ipw.BoundedFloatText(min=0, max=1, step=0.01)
        _check = ipw.Checkbox(indent=False, layout=ipw.Layout(width='30px'))
        _item.index = self.no_of_elements - 1
        _check.index = self.no_of_elements - 1
        _item.observe(self.input_item_eventhandler, "value")
        _check.observe(self.check_eventhandler, "value")
        self.input_items.append(_item)
        self.checks.append(_check)
        
        # Keep the order of description same with atomic number
        self.elements = self.get_symbols()
        for i in range(self.no_of_elements):
            self.input_items[i].description = self.elements[i]

        self.box_abundance_editor.children = [ipw.VBox(self.input_items), ipw.VBox(self.checks)]
        
        # Add new trace to plot.
        self.fig.add_scatter(x=self.velocity[1:], # convert to km/s
                        y=[0]*self.no_of_shells,
                        mode="lines+markers",
                        name=element_symbol_string,
                    )
        # Sort the legend in atomic order.
        fig_data_lst = list(self.fig.data)
        fig_data_lst.insert(np.argwhere(self.elements == element_symbol_string)[0][0]+1, self.fig.data[-1])
        self.fig.data = fig_data_lst[:-1]
        
        self.read_abundance()
            
        # Clear symbol input box.
        self.input_symb.value = ''
        
    



        

    @out_plot.capture(clear_output=True)
    def generate_abundance_density_plot(self):
        self.fig = go.FigureWidget()
        title = "Abundance/Density vs Velocity"
        data = self.abundance.T
        
        self.fig.add_trace(
            go.Scatter(
                x=self.velocity[1:],
                y=self.density[1:],
                mode="lines+markers",
                name="<b>Density</b>",
                yaxis="y2",
                line=dict(color="black"),
                marker_symbol="square",
            ),
        )
        
        for i in range(self.no_of_elements):
            self.fig.add_trace(
                go.Scatter(
                    x=self.velocity[1:],
                    y=data.iloc[:,i],
                    mode="lines+markers",
                    name=self.elements[i],
                ),
            )
            
        self.fig.update_layout(
            xaxis=dict(
                title="Velocity (km/s)",
                tickformat = "f"
            ),
            yaxis=dict(title="Fractional Abundance", 
                    exponentformat="e",
                    range=[0, 1]
                    ),
            yaxis2=dict(title="Density", 
                        exponentformat="e",                    
                        overlaying="y",
                        side="right"),
            shapes=[         # Line Diagonal
                    dict(
                        type="line",
                        yref="y",
                        y0=0,
                        y1=1,
                        xref="x",
                        x0=self.velocity[1].value,
                        x1=self.velocity[1].value,
                        line=dict(
                            color="Red",
                            dash="dot",
                        )
                    )],
            height=500,
            title=title,
            hovermode="closest",
            legend=dict(
                x=1.15,
            )
        )
            


    def display(self):
        self.box_abundance_editor = ipw.HBox([ipw.VBox(self.input_items), ipw.VBox(self.checks)])


        help_note = ipw.HTML(value="<p style=\"text-indent: 40px\"><b>Click the checkbox</b> to </p> <p style=\"text-indent: 40px\"> 1) lock the abundance you don't want to normalize </p> <p style=\"text-indent: 40px\"> 2) apply the abundance to other shells.</p>",
                        indent=True
                        )



        norm_box = ipw.HBox([self.btn_norm, self.norm_warning])

    @classmethod
    def from_csvy(cls, fpath):
        csvy_model_config, csvy_model_data = load_csvy(fpath)
        base_dir = os.path.abspath(os.path.dirname(__file__))
        schema_dir = os.path.join(base_dir, "../..", "io", "schemas")
        csvy_schema_file = os.path.join(schema_dir, "csvy_model.yml")
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_file)
        )

        if hasattr(csvy_model_config, "velocity"):
            velocity = quantity_linspace(csvy_model_config.velocity.start,
                                         csvy_model_config.velocity.stop,
                                         csvy_model_config.velocity.num + 1).cgs
        else:
            velocity_field_index = [field["name"] for field in csvy_model_config.datatype.fields].index("velocity")
            velocity_unit = u.Unit(csvy_model_config.datatype.fields[velocity_field_index]["unit"])
            velocity = csvy_model_data["velocity"].values * velocity_unit

        no_of_shells = len(velocity) - 1

        if hasattr(csvy_model_config, "density"):
            adjusted_velocity = velocity.insert(0, 0)
            v_middle = (adjusted_velocity[1:] * 0.5 +
                        adjusted_velocity[:-1] * 0.5)
            no_of_shells = len(adjusted_velocity) - 1

            d_conf = csvy_model_config.density
            density_type = d_conf.type
            if density_type == "branch85_w7":
                density_0 = calculate_power_law_density(v_middle, d_conf.w7_v_0,
                                                        d_conf.w7_rho_0, -7)
            elif density_type == "uniform":
                density_0 = (d_conf.value.to("g cm^-3") *
                             np.ones(no_of_shells))
            elif density_type == "power_law":
                density_0 = calculate_power_law_density(v_middle, d_conf.v_0,
                                                        d_conf.rho_0,
                                                        d_conf.exponent)
            elif density_type == "exponential":
                density_0 = calculate_exponential_density(v_middle, d_conf.v_0,
                                                          d_conf.rho_0)
            else:
                raise ValueError(f"Unrecognized density type "
                                 f"{d_conf.type}")
        else:
            density_field_index = [field["name"] for field in csvy_model_config.datatype.fields].index("density")
            density_unit = u.Unit(csvy_model_config.datatype.fields[density_field_index]["unit"])
            density_0 = csvy_model_data["density"].values * density_unit

        if hasattr(csvy_model_config, "abundance"):
            abundances_section = csvy_model_config.abundance
            abundance, isotope_abundance = read_uniform_abundances(
                abundances_section, no_of_shells
            )
        else:
            _, abundance, isotope_abundance = parse_csv_abundances(
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
            density=density_0,
            abundance=abundance, 
            velocity=velocity
            )

    @classmethod
    def from_yml(cls, fpath):
        config = Configuration.from_yaml(fpath)

        if hasattr(config, "csvy_model"):
            model = Radial1DModel.from_csvy(config)
        else:
            model = Radial1DModel.from_config(config)

        velocity = model.velocity
        density = model.homologous_density.density_0
        abundance = model.abundance
        isotope_abundance = model.isotope_abundance

        # Combine elements and isotopes to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        return cls(
            density=density,
            abundance=abundance, 
            velocity=velocity
            )

    
    @classmethod
    def from_hdf(cls, fpath):
        with pd.HDFStore(fpath, "r") as hdf:
            abundance = hdf["/simulation/plasma/abundance"]
            _density = hdf["/simulation/model/homologous_density/density_0"]
            v_inner = hdf["/simulation/model/v_inner"]
            v_outer = hdf["/simulation/model/v_outer"]
        
        density = np.array(_density) * u.g / (u.cm)**3
        velocity = np.append(v_inner, v_outer[len(v_outer)-1]) * u.cm / u.s

        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        
        return cls(density=density, abundance=abundance, velocity=velocity)


    @classmethod
    def from_simulation(cls, sim):
        abundance = sim.model.raw_abundance.copy()
        isotope_abundance = sim.model.raw_isotope_abundance.copy()

        # integrate element and isotope to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        velocity = sim.model.velocity
        density = sim.model.homologous_density.density_0

        return cls(density=density, abundance=abundance, velocity=velocity)



