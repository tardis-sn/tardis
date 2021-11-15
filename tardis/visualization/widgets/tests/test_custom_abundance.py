"""Tests for custom abundance widget."""
import os
import pytest
import tardis
import numpy as np
import numpy.testing as npt
import pandas.testing as pdt
import plotly.graph_objects as go
from astropy import units as u

from tardis.visualization.widgets.custom_abundance import(
    CustomAbundanceWidgetData,
    CustomYAML,
    CustomAbundanceWidget,
    DensityEditor
)

@pytest.fixture(scope="module")
def yml_data():
    yml_path = os.path.join(tardis.__path__[0], "io", "tests", "data", "tardis_configv1_verysimple.yml")
    return CustomAbundanceWidgetData.from_yml(yml_path)

@pytest.fixture(scope="module")
def csvy_data():
    csvy_path = os.path.join(tardis.__path__[0], "io", "tests", "data", "csvy_full.csvy")
    return CustomAbundanceWidgetData.from_csvy(csvy_path)

@pytest.fixture(scope="module")
def hdf_data(hdf_file_path, simulation_verysimple):
    simulation_verysimple.to_hdf(
        hdf_file_path, overwrite=True
    )  # save sim at hdf_file_path
    return CustomAbundanceWidgetData.from_hdf(hdf_file_path)

@pytest.fixture(scope="module")
def sim_data(simulation_verysimple):
    return CustomAbundanceWidgetData.from_simulation(simulation_verysimple)

# @pytest.fixture(scope="module", params=[yml_data, csvy_data, hdf_data, sim_data])
@pytest.fixture(scope="module")
def caw():
    return CustomAbundanceWidget(yml_data)

class TestCustomAbundanceWidgetData:
    def test_get_symbols(self):
        symbols = caw.get_symbols()
        npt.assert_array_equal(symbols, ["O", "Mg", "Si", "S", "Ar", "ca"])

class TestCustomAbundanceWidget:
    def test_update_input_item_value(self):
        caw.update_input_item_value(0, 0.33333)
        assert caw.input_items[0].value == 0.333

    def test_read_abundance(self):
        caw.data.abundance[0] = 0.2
        caw.read_abundance()
        for i in range(caw.no_of_elements):
            assert caw.input_items[i].value == 0.2

    def test_update_abundance_plot():
        caw.data.abundance.iloc[0, :] = 0.2
        caw.update_abundance_plot(0)

        npt.assert_array_equal(caw.fig.data[2].y, np.array([0.2]*(caw.no_of_shells+1)))

    def test_bound_locked_sum_to_1(self):
        """Trigger checkbox eventhandler and input_item eventhandler 
        to test `bound_locked_sum_to_1()` function
        """
        # bound checked input to 1
        caw.checks[0].value = True
        caw.checks[1].value = True
        caw.input_items[0].value = 0.5
        caw.input_items[1].value = 0.6
        assert caw.input_items[1].value == 0.5

        # bound to 1 when input is checked
        caw.checks[2].value = True
        assert caw.checks[2].value == 0

    @pytest.mark.parametrize(("v0", "v1", "expected"), [(11000, 11450, "hidden"), (11100, 11200, "hidden"), (11000, 11451, "visible")])
    def test_overwrite_existing_shells(self, v0, v1, expected):
        caw.input_v_start.value = v0
        caw.input_v_end.value = v1

        assert caw.overwrite_warning.layout.visibility == expected
        
    def test_update_bar_diagonal():
        pass

    def test_on_btn_norm():
        pass

    def test_on_btn_add_element():
        pass

    def test_on_btn_add_shell():
        pass


class TestCustomYAML:
    def test_create_fields_dict(self):
        custom_yaml = CustomYAML("test", 0, 0, 0, 0)
        custom_yaml.create_fields_dict(["H", "He"])
        datatype_dict = {
            "fields":[
                {
                    "name": "velocity",
                    "unit": "km/s"
                },
                {
                    "name": "density",
                    "unit": "g/cm^3"
                },
                {
                    "name": "H",
                    "desc": "fractional H abundance"
                },
                {
                    "name": "He",
                    "desc": "fractional He abundance"
                }
            ]
        }

        assert custom_yaml.datatype == datatype_dict

class TestDensityEditor:
    pass