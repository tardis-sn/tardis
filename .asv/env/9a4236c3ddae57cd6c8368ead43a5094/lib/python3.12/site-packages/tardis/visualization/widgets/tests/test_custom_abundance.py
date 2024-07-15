"""Tests for custom abundance widget."""
from pathlib import Path
import pytest
import tardis
import numpy as np
import numpy.testing as npt

from tardis.tests.test_util import monkeysession
from tardis.visualization.widgets.custom_abundance import (
    CustomAbundanceWidgetData,
    CustomYAML,
    CustomAbundanceWidget,
)


@pytest.fixture(scope="module")
def yml_data(example_configuration_dir: Path, atomic_dataset):
    """Fixture to contain a CustomAbundanceWidgetData
    instance generated from a YAML file tardis_configv1_verysimple.yml.

    Returns
    -------
    CustomAbundanceWidgetData
        CustomAbundanceWidgetData generated from a YAML
    """
    yml_path = example_configuration_dir / "tardis_configv1_verysimple.yml"

    return CustomAbundanceWidgetData.from_yml(
        yml_path, atom_data=atomic_dataset
    )


@pytest.fixture(scope="module")
def caw(yml_data, monkeysession):
    """Fixture to contain a CustomAbundanceWidget
    instance generated from a YAML file tardis_configv1_verysimple.yml.

    Returns
    -------
    CustomAbundanceWidget
        CustomAbundanceWidget generated from a YAML
    """
    caw = CustomAbundanceWidget(yml_data)
    monkeysession.setattr(
        "tardis.visualization.widgets.custom_abundance.is_notebook",
        lambda: True,
    )
    caw.display()
    return caw


class TestCustomAbundanceWidgetData:
    def test_get_symbols(self, yml_data):
        """Tests the atomic symbols for the YAML CustomAbundanceWidgetData"""
        symbols = yml_data.get_symbols()
        npt.assert_array_equal(symbols, ["O", "Mg", "Si", "S", "Ar", "Ca"])


class TestCustomAbundanceWidget:
    def test_update_input_item_value(self, caw):
        """Tests updating an input item value

        Parameters
        ----------
        caw : CustomAbundanceWidget
        """
        caw.update_input_item_value(0, 0.33333)
        assert caw.input_items[0].value == 0.333

    def test_read_abundance(self, caw):
        """Tests reading an abundance

        Parameters
        ----------
        caw : CustomAbundanceWidget
        """
        caw.data.abundance[0] = 0.2
        caw.read_abundance()
        for i in range(caw.no_of_elements):
            assert caw.input_items[i].value == 0.2

    def test_update_abundance_plot(self, caw):
        """Tests plotting an abundance array

        Parameters
        ----------
        caw : CustomAbundanceWidget
        """
        caw.data.abundance.iloc[0, :] = 0.2
        caw.update_abundance_plot(0)

        npt.assert_array_equal(
            caw.fig.data[2].y, np.array([0.2] * (caw.no_of_shells + 1))
        )

    def test_bound_locked_sum_to_1(self, caw):
        """Trigger checkbox eventhandler and input_item eventhandler
        to test `bound_locked_sum_to_1()` function.
        """
        # bound checked input to 1
        caw.checks[0].value = True
        caw.checks[1].value = True
        caw.input_items[0].value = 0.5
        caw.input_items[1].value = 0.6
        assert caw.input_items[1].value == 0.5

        # bound to 1 when input is checked
        caw.checks[2].value = True
        assert caw.input_items[2].value == 0

    @pytest.mark.parametrize(
        "v0, v1, expected",
        [
            (11000, 11450, "hidden"),
            (11100, 11200, "hidden"),
            (11000, 11451, "visible"),
        ],
    )
    def test_overwrite_existing_shells(self, caw, v0, v1, expected):
        """Trigger velocity input box handler to test whether overwriting
        existing shell.
        """
        caw.input_v_start.value = v0
        caw.input_v_end.value = v1

        assert caw.overwrite_warning.layout.visibility == expected

    @pytest.mark.skip(
        reason="Problem with cutting model in restructure. See TARDIS issue 2390"
    )
    @pytest.mark.parametrize(
        "multishell_edit, expected_x, expected_y, expected_width",
        [
            (False, [19775], [1], [450]),
            (True, (11225, 15500), (1, 1), (450, 9000)),
        ],
    )
    def test_update_bar_diagonal(
        self, caw, multishell_edit, expected_x, expected_y, expected_width
    ):
        """Tests updating the bar figure

        Parameters
        ----------
        caw : CustomAbundanceWidget
        multishell_edit : bool
        expected_x : list
        expected_y : list
        expected_width : list
        """
        if multishell_edit:
            caw.irs_shell_range.disabled = False  # update_bar_diagonal() will be called when status of irs_shell_range is changed
            caw.irs_shell_range.value = (1, 20)

            assert caw.shell_no == caw.irs_shell_range.value[0]
            assert caw.btn_next.disabled == True
            assert caw.btn_prev.disabled == True
        else:
            caw.shell_no = 20
            caw.update_bar_diagonal()

        npt.assert_array_almost_equal(caw.fig.data[0].x, expected_x)
        npt.assert_array_almost_equal(caw.fig.data[0].width, expected_width)
        npt.assert_array_almost_equal(caw.fig.data[0].y, expected_y)

    @pytest.mark.parametrize(
        "multishell_edit, inputs, locks, expected",
        [
            (False, [0, 0, 0, 0, 0, 0, 0], [False] * 6, [0, 0, 0, 0, 0, 0, 0]),
            (
                False,
                [0.1, 0.2, 0, 0, 0, 0, 0],
                [True] + [False] * 5,
                [0.1, 0.9, 0, 0, 0, 0, 0],
            ),
            (
                False,
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                [False] * 6,
                [0.0476, 0.0952, 0.143, 0.19, 0.238, 0.286],
            ),
            (
                False,
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                [True] * 2 + [False] * 4,
                [0.1, 0.2, 0.117, 0.156, 0.194, 0.233],
            ),
            (
                True,
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                [True] * 2 + [False] * 4,
                [0.1, 0.2, 0.117, 0.156, 0.194, 0.233],
            ),
        ],
    )
    def test_on_btn_norm(self, caw, multishell_edit, inputs, locks, expected):
        """Tests the normalisation button

        Parameters
        ----------
        caw : CustomAbundanceWidget
        multishell_edit : bool
        inputs : list
        locks : list
        expected : list
        """
        if multishell_edit:
            caw.rbs_multi_apply.index = 0
            for i, item in enumerate(caw.input_items):
                item.value = inputs[i]
                caw.checks[i].value = locks[i]

            caw.on_btn_norm(None)

            for i, item in enumerate(caw.input_items):
                assert item.value == expected[i]

            start_no = caw.irs_shell_range.value[0]
            end_no = caw.irs_shell_range.value[1]

            for i, v in enumerate(expected):
                line = caw.fig.data[2 + i].y[start_no - 1 : end_no]
                unique_v = set(line)
                assert len(unique_v) == 1
                unique_v = float("{:.3g}".format(list(unique_v)[0]))
                assert unique_v == v
        else:
            for i, item in enumerate(caw.input_items):
                item.value = inputs[i]
                caw.checks[i].value = locks[i]

            caw.on_btn_norm(None)

            for i, item in enumerate(caw.input_items):
                assert item.value == expected[i]


class TestCustomYAML:
    def test_create_fields_dict(self):
        """Test creating fields in the YAML"""
        custom_yaml = CustomYAML("test", 0, 0, 0, 0)
        custom_yaml.create_fields_dict(["H", "He"])
        datatype_dict = {
            "fields": [
                {"name": "velocity", "unit": "km/s"},
                {"name": "density", "unit": "g/cm^3"},
                {"name": "H", "desc": "fractional H abundance"},
                {"name": "He", "desc": "fractional He abundance"},
            ]
        }

        assert custom_yaml.datatype == datatype_dict
