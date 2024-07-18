from tardis.io.configuration.config_reader import ConfigurationNameSpace
import pytest
from astropy import units as u
import os

from numpy.testing import assert_almost_equal

simple_config_dict = {
    "a": {"b": {"param1": 1, "param2": [0, 1, 2 * u.km], "param3": 4.0 * u.km}}
}


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, "data", filename)


@pytest.fixture(scope="function")
def config_ns():
    """Return example instance of `ConfigurationNameSpace` class."""
    config_name_space = ConfigurationNameSpace(simple_config_dict)
    return config_name_space


def test_simple_configuration_namespace(config_ns):
    assert config_ns.a.b.param1 == 1
    config_ns.a.b.param1 = 2
    assert config_ns["a"]["b"]["param1"] == 2

    config_ns["a"]["b"]["param1"] = 3
    assert config_ns.a.b.param1 == 3


def test_quantity_configuration_namespace(config_ns):
    config_ns.a.b.param3 = 3
    assert_almost_equal(config_ns["a"]["b"]["param3"].to(u.km).value, 3)

    config_ns.a.b.param3 = 5000 * u.m
    assert_almost_equal(config_ns["a"]["b"]["param3"].to(u.km).value, 5)


def test_access_with_config_item_string(config_ns):
    assert config_ns.get_config_item("a.b.param1") == 1

    config_ns.set_config_item("a.b.param1", 2)
    assert config_ns.a.b.param1 == 2


def test_set_with_config_item_string_quantity(config_ns):
    config_ns.set_config_item("a.b.param3", 2)
    assert_almost_equal(config_ns.a.b.param3.to(u.km).value, 2)


def test_get_with_config_item_string_item_access(config_ns):
    item = config_ns.get_config_item("a.b.param2.item0")
    assert item == 0
    item = config_ns.get_config_item("a.b.param2.item1")
    assert item == 1


def test_set_with_config_item_string_item_access(config_ns):
    config_ns.set_config_item("a.b.param2.item0", 2)

    item = config_ns.get_config_item("a.b.param2.item0")

    assert item == 2


def test_set_with_config_item_string_item_access_quantity(config_ns):
    config_ns.set_config_item("a.b.param2.item2", 7)

    item = config_ns.get_config_item("a.b.param2.item2")

    assert_almost_equal(item.to(u.km).value, 7)


def test_config_namespace_copy(config_ns):
    config_ns2 = config_ns.deepcopy()
    config_ns2.a.b.param1 = 2
    assert config_ns2.a.b.param1 != config_ns.a.b.param1


def test_config_namespace_quantity_set():
    data_path("paper1_tardis_configv1.yml")
