import os
from glob import glob

from astropy import units as u
import pytest

from tardis.io.config_validator import DefaultParser, Config, ConfigValueError


existing_configs = glob(os.path.join('docs', 'examples', '*.yml'))
existing_configs += glob(os.path.join('tardis', 'io', 'tests', 'data', '*.yml'))
config_definition = os.path.join('tardis', 'data', 'tardis_config_definition.yml')

test_config_definition = os.path.join('tardis', 'io', 'tests', 'data', 'conf_def.yml')
test_config = os.path.join('tardis', 'io', 'tests', 'data', 'conf_tes.yml')
existing_configs.remove(test_config_definition)
existing_configs.remove(test_config)


@pytest.mark.parametrize("config_filename", existing_configs)
def test_configread(config_filename):
    config = Config.from_yaml(config_filename, config_definition)


def test_configread_test_config():
    config = Config.from_yaml(test_config, test_config_definition)


def default_parser_helper(test_dic, default, wdefault, value, wvalue, container, mandatory, return_default=None,
                          return_value=None, value_as_string=True):
    test_ob = DefaultParser(test_dic)

    if return_value is None:
        return_value = value

    if return_default is None:
        return_default = default

    if not default == None:
        dhelper = True
    else:
        dhelper = False

    assert test_ob.has_default == dhelper
    assert test_ob.get_default() == return_default
    assert test_ob.is_leaf
    assert test_ob.is_container() == container
    assert test_ob.is_mandatory == mandatory

    #set good default
    test_ob.set_default(default)

    assert test_ob.get_default() == return_default
    #set bad default
    if wdefault is not None:
        with pytest.raises(ValueError):
            test_ob.set_default(wdefault)

    assert test_ob.get_value() == return_default

    #set good value
    test_ob.set_config_value(value)

    assert test_ob.get_value() == return_value

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

    default = "99.99 cm"
    return_default = 99.99 * u.cm
    wdefault = "kl"
    value = "11.12 m"
    return_value = 11.12 * u.m
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)




def test_default_parser_quantity_range():
    example_dic = {'default': ['1 cm', '5 cm'],
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'quantity_range'}

    default = ['1.0 cm', '5 cm']
    return_default = [1.0 * u.cm, 5 * u.cm]
    wdefault = "kl"
    value = ['10 m', '50 cm']
    return_value = [10 * u.m, 50 * u.cm]
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)


def test_default_parser_quantity_range_old():
    example_dic = {'default': {'start': '1 cm', 'end': '5 cm'},
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'quantity_range'}

    default = ['1.0 cm', '5 cm']
    return_default = [1.0 * u.cm, 5 * u.cm]
    wdefault = "kl"
    value = ['10 m', '50 cm']
    return_value = [10 * u.m, 50 * u.cm]
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)


def test_default_parser_quantity_range_sampeled():
    example_dic = {'default': ['1 cm', '5 cm', 10],
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'quantity_range_sampled'}

    default = ['1.0 cm', '5 cm', 10]
    return_default = [1.0 * u.cm, 5 * u.cm, 10]
    wdefault = "kl"
    value = ['10 m', '50 cm', 10]
    return_value = [10 * u.m, 50 * u.cm, 10]
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)


def test_default_parser_quantity_range_sampeled_old():
    example_dic = {'default': {'start': '1 cm', 'stop': '5 cm', 'num': 10},
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'quantity_range_sampled'}

    default = ['1.0 cm', '5 cm', 10]
    return_default = [1.0 * u.cm, 5 * u.cm, 10]
    wdefault = "kl"
    value = ['10 m', '50 cm', 10]
    return_value = [10 * u.m, 50 * u.cm, 10]
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)


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


def test_default_parser_range_old():
    example_dic = {'default': {'start': 0, 'stop': 10},
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


def test_default_parser_range_sampled():
    example_dic = {'default': {'start': 0, 'stop': 10, 'num': 1},
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


def test_property_type_abundances():
    example_dic = {'default': {'He': 0.4, 'Mg': 0.1, 'Pb': 0.5},
                   'help': 'quantity for testing',
                   'mandatory': True,
                   'property_type': 'abundance_set'}

    default = {'He': 0.4, 'Mg': 0.1, 'Pb': 0.5}
    return_default = {'he': 0.4, 'mg': 0.1, 'pb': 0.5}
    wdefault = "kl"
    value = {'He': 0.4, 'Mg': 0.6}
    return_value = {'he': 0.4, 'mg': 0.6}
    wvalue = "yy"
    container = False
    mandatory = True

    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory,
                               return_default=return_default, return_value=return_value)






