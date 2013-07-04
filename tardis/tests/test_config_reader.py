# tests for the config reader module

from tardis import config_reader
from astropy import units as u
import pytest



def test_quantity_parser_normal():
    q1 = config_reader.parse2quantity('5 km/s')
    assert q1.value == 5.
    assert q1.unit == u.Unit('km/s')

def test_quantity_parser_malformed_quantity1():
    with pytest.raises(config_reader.TardisConfigError):
        q1 = config_reader.parse2quantity('abcd')

def test_quantity_parser_malformed_quantity2():
    with pytest.raises(config_reader.TardisMalformedQuantityError):
        q1 = config_reader.parse2quantity('5 abcd')


def test_config_namespace_attribute_test():
    namespace = config_reader.TardisConfigurationNameSpace({'param1':1})
    assert namespace.param1 == 1

def test_config_namespace_attribute_test():
    namespace = config_reader.TardisConfigurationNameSpace({'param1':1})
    with pytest.raises(AttributeError):
        assert namespace.param2 == 1
