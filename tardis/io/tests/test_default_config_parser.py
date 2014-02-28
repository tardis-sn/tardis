
import tardis
from tardis.io.default_config_parser import *
import os
import yaml





#test the whole thing
def test_default_config_parser():
    test_config_fname = os.path.join(tardis.__path__[0], 'data', 'conf_tes.yml')
    test_default_fname = os.path.join(tardis.__path__[0], 'data', 'conf_def.yml')

    with open(test_config_fname) as f:
        test_config = yaml.SafeLoader(f)

    with open(test_default_fname) as f:
        test_default = yaml.SafeLoader(f)
    

    test_conf_ob = Config(test_default, test_config)


def default_parser_helper(test_dic, default, wdefault, value, wvalue, container, mandatory):
    test_ob = DefaultParser(test_dic)

    if not default == None:
        dhelper  = True
    else:
        dhelper = False
        
    if not test_ob.has_default() == dhelper:
        return 1
    print(test_ob.get_default())
    print(default)
    if not test_ob.get_default() == default:
        return 2
    if not test_ob.is_leaf:
        return 3
    if not test_ob.is_container() == container:
        return 4
    if not test_ob.is_mandatory() == mandatory:
        return 5
    
    #set good default
    try:
        test_ob.set_default(default)
    except:
        return 6
    
    if not test_ob.get_default() == default:
        return 7
    #set bad default
    if not wdefault == None:
        try:
            test_ob.set_default(wdefault)
            return 8
        except ValueError:
            pass
    
    if not test_ob.get_value() == default:
        return 9
    
    #set good value
    try:
        test_ob.set_config_value(value)
    except:
        return 10
    
    if not test_ob.get_value() == value:
        return 11
    
    #set bad value
    if not wvalue == None:
        try:
            test_ob.set_config_value(wvalue)
            if test_ob.get_value() == wvalue:
                return 12
            elif test_ob.get_value() == value:
                return 12
        except:
            pass
    
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
    if not  ex == 0:
        print ex
        raise AssertionError
    
    
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
    if not  ex == 0:
        print ex
        raise AssertionError
    
def test_default_parser_quantity():
    example_dic = {'default': '99.99 cm',
    'help': 'quantity for testing',
    'mandatory': True,
    'property_type': 'quantity'}
    
    default = "99.99 cm"
    wdefault = "kl"
    value = "11.12 m"
    wvalue = "yy"
    container = False
    mandatory = True
    
    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)
    if not  ex == 0:
        print ex
        raise AssertionError
    
    
def test_default_parser_range():
    example_dic = {'default': [0, 10],
    'help': 'range for testing',
    'mandatory': False,
    'property_type': 'range'}
    
    default = [0,10]
    wdefault = 1
    value = [7,8]
    wvalue = 2
    container = False
    mandatory = False
    
    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)
    if not  ex == 0:
        print ex
        raise AssertionError

def test_default_parser_range_sampled():
    example_dic = {'default': [0, 10, 1],
    'help': 'range for testing',
    'mandatory': False,
    'property_type': 'range_sampled'}
    
    default = [0,10,1]
    wdefault = [1,3]
    value = [1,5,1]
    wvalue = [1,1]
    container = False
    mandatory = False
    
    ex = default_parser_helper(example_dic, default, wdefault, value, wvalue, container, mandatory)
    if not  ex == 0:
        print ex
        raise AssertionError
    
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
    if not  ex == 0:
        print ex
        raise AssertionError
 
