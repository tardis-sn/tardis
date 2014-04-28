import os
from glob import glob

from astropy import units as u
import pytest

from tardis.io.default_config_parser import DefaultParser, Config, ConfigValueError

existing_configs = glob(os.path.join('tardis', 'docs', 'examples', '*.yml'))
config_definition = os.path.join('tardis', 'data', 'tardis_default_config_definition.yml')


@pytest.mark.parametrize(("config_filename",), existing_configs)
def test_configread(config_filename):
    config = Config.from_yaml(config_filename, config_definition)


def default_parser_helper(test_dic, default, wdefault, value, wvalue, container, mandatory):
    test_ob = DefaultParser(test_dic)

    if not default == None:
        dhelper = True
    else:
        dhelper = False

    assert test_ob.has_default == dhelper
    assert test_ob.get_default() == default
    assert test_ob.is_leaf
    assert test_ob.is_container() == container
    assert test_ob.is_mandatory == mandatory

    #set good default
    test_ob.set_default(str(default))

    assert test_ob.get_default() == default
    #set bad default
    if wdefault is not None:
        with pytest.raises(ValueError):
            test_ob.set_default(wdefault)

    assert test_ob.get_value() == default

    #set good value
    test_ob.set_config_value(str(value))

    assert test_ob.get_value() == value

    #set bad value
    if wvalue is not None:
        with pytest.raises(ConfigValueError):
            test_ob.set_config_value(wvalue)
            test_ob.get_value()

    return 0


def test_default_parser_float():
    example_dic = {'default': 99.99,
                   'help': 'float value for testing',
                   'mandatory': True,
                   'property_type': 'float'}
    default = 99.99
    wdefault = "xx"
    value = 11.12
    wvalue = "yy"
    container = False
    mandatory = True
    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)


def test_default_parser_integer():
    example_dic = {'default': 99,
                   'help': 'integer value for testing',
                   'mandatory': True,
                   'property_type': 'int'}

    default = 99
    wdefault = 9.15
    value = 11
    wvalue = 9.22
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)


def test_default_parser_quantity():
    example_dic = {'default': '99.99 cm',
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'quantity'}

    default = 99.99 * u.cm
    wdefault = "kl"
    value = 11.12 * u.m
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)

# def test_default_parser_quantity_range():
#     example_dic = {'default': ['1 cm', '5 cm'],
#                    'help': 'quantity for testing',
#                    'mandatory': True,
#                    'property_type': 'quantity_range'}
#
#     default = [1.0 * u.cm, 5 * u.cm]
#     wdefault = "kl"
#     value = [10 * u.m, 50 * u.cm]
#     wvalue = "yy"
#     container = False
#     mandatory = True
#
#     ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)

def test_default_parser_range():
    example_dic = {'default': [0, 10],
                   'help': 'range for testing',
                   'mandatory': False,
                   'property_type': 'range'}

    default = [0, 10]
    wdefault = 1
    value = [7, 8]
    wvalue = 2
    container = False
    mandatory = False

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)


def test_default_parser_range_sampled():
    example_dic = {'default': [0, 10, 1],
                   'help': 'range for testing',
                   'mandatory': False,
                   'property_type': 'range_sampled'}

    default = [0, 10, 1]
    wdefault = [1, 3]
    value = [1, 5, 1]
    wvalue = [1, 1]
    container = False
    mandatory = False

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)


def test_default_parser_string():
    example_dic = {'default': 'DEFAULT',
                   'help': 'string for testing',
                   'mandatory': True,
                   'property_type': 'string'}

    default = "DEFAULT"
    wdefault = None
    value = "blub"
    wvalue = None
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)
 
