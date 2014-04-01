# coding=utf-8

import re
import logging
from astropy import units
from tardis.util import parse_quantity
from astropy.units.core import UnitsException
import yaml

logger = logging.getLogger(__name__)



class Error(Exception):
    """Base class for exceptions in the config parser."""
    pass


class ConfigTypeError(Error):
    """
    Exception raised if the type of the configured value mismatches the type
    specified in the default configuration.
    """

    def __init__(self, value, expected_type, help):
        self.value = value
        self.expected_type = expected_type
        self.help = help

    def __str__(self):
        return "Expected type %s but found %s.\nHelp:%s " % \
        (repr(self.expected_type), repr(type(self.value)), help)


class ConfigError(Error):
    """
    Exception raised if something is wrong in the default configuration.
    """

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return "Error in the configuration at %s " % ("->".join(self.path))


class DefaultConfigError(ConfigError):
    """
    Exception raised if something is wrong in the default configuration.
    """

    def __str__(self):
        return "Error in the default configuration at %s " % \
        ("->".join(self.path))
    
    
class PropertyType(object):
    def __init__(self):
        self._default = None
        self._allowed_value = None
        self._allowed_type = None
        self._help = None
        self._mandatory = False
        self._lower = None
        self._upper =None
        pass
    
    @property
    def default(self):
        return self._default
    
    @default.setter
    def default(self, value):
        self._default = self.to_type(value)
        
    @property
    def allowed_value(self):
        return self._allowed_value
    
    @allowed_value.setter
    def allowed_value(self, value):
        if isinstance(value, basestring):
            self._allowed_value = set(self.__list_dtype(value.split()))
        elif isinstance(value, list) or isinstance(value, set):
            self._allowed_value = set(self.__list_dtype(value))
        elif isinstance(value, float) or isinstance(value, int):
            self._allowed_type = set([value])
        else:
            raise ValueError("Can not set allowed value.")
            
            
    @property
    def allowed_type(self):
        return self._allowed_value
    
    @allowed_type.setter
    def allowed_type(self, value):
        self._allowed_type = value
        if '_parse_allowed_type' in (set(dir(self.__class__)) - set(dir(PropertyType))) and value != None:
            self._lower, self._upper = self._parse_allowed_type(value)
        
    @property
    def help(self):
        return self._help
    
    @help.setter
    def help(self, value):
        self._help = value
        
    @property
    def mandatory(self):
        return self._mandatory
    
    @mandatory.setter
    def mandatory(self, value):
        self._mandatory = value
    
    def check_type(self, value):
        return True
        
    def to_type(self, value):
        return value
    
    def __list_dtype(self, mixed_list):
        try:
            tmp = [str(a) for a in mixed_list]
        except ValueError:
            try:
                tmp = [float(a) for a in mixed_list]
            except ValueError:
                try:
                    tmp = [int(a) for a in mixed_list]
                except:
                    raise ValueError("Forbidden type in allowed_type")
        return tmp
        
    
    def _check_allowed_value(self, _value):
        """
        Returns True if the value is allowed or no allowed value is given.
        """
        if self._allowed_value != None:
            atype = type(iter(self.allowed_value).next())
            value = atype(_value)
            if value in self.allowed_value:
                return True
            else:
                return False
        else:
            return True
        
 
class PropertyTypeContainer(PropertyType):
    
    def check_type(self):
        pass
    
   
class PropertyTypeInt(PropertyType):
    
    def check_type(self, value):
        try:
            int(value)
            if float.is_integer(float(value)):
                if self._check_allowed_value(value) and self._check_allowed_type(value):
                    return True
            else:
                return False
        except ValueError:
            return False
    
    def to_type(self, value):
        return int(value)
 #ToDo: use this if allowed type is specified    
    def _parse_allowed_type(self, allowed_type):
        string = allowed_type.strip()
        upper = None
        lower = None
        if string.find("<") or string.find(">"):
            #like x < a
            match = re.compile('[<][\s]*[0-9.+^*eE]*$').findall(string)
            if match:
                value = re.compile('[0-9.+^*eE]+').findall(string)[0]
                upper = float(value)
                #like a > x"
            match = re.compile('^[\s0-9.+^*eE]*[\s]*[<]$').findall(string)
            if match:
                value = re.compile('[0-9.+^*eE]+').findall(string)[0]
                upper = float(value)
                #like x > a
            match = re.compile('[>][\s]*[0-9.+^*eE]*$').findall(string)
            if match:
                value = re.compile('[0-9.+^*eE]+').findall(string)[0]
                lower = float(value)
                #like a < x
            match = re.compile('^[\s0-9.+^*eE]*[\s]*[<]$').findall(string)
            if match:
                value = re.compile('[0-9.+^*eE]+').findall(string)[0]
                lower = float(value)
        return lower, upper
    
    def __check_type(self, value, lower_lim, upper_lim):
        upper, lower = True, True
        if upper_lim != None:
            upper = value < upper_lim
        if lower_lim != None:
            lower = value > lower_lim
        return upper and lower
    
    
    def _check_allowed_type(self, value):
        if self._allowed_type != None:
            if  self.__check_type(value, self._lower, self._upper):
                return True
            else:
                return False
        else:
            return True
    
    def _is_valid(self, value):
        if not self.check_type(value):
            return False
        if self.allowed_value != None:
            return False
        if not self.__check_type(value, self._lower, self._upper):
            return  False
        return True
            
            
        
    
class PropertyTypeFloat(PropertyTypeInt):
    
    def check_type(self, value):
        try:
            float(value)
            if self._check_allowed_value(value) and self._check_allowed_type(value):
               return True
        except ValueError:
            return False
        
    def to_type(self, value):
        return float(value)
    
    
class PropertyTypeQuantity(PropertyType):
    
    def check_type(self, value):
        try:
            quantity_value, quantity_unit = value.strip().split()
            float(quantity_value)
            units.Unit(quantity_unit)
            return True
        except ValueError:
            return False
    
    def to_type(self, value):
        quantity_value, quantity_unit = value.strip().split()
        float(quantity_value)
        units.Unit(quantity_unit)
        return (quantity_value, quantity_unit)
    
class PropertyTypeQuantityRange(PropertyTypeQuantity):
    
    def _to_units(self, los):
        if len(los) > 2:
            loq = [(lambda x: (units.Quantity(float(x[0]),x[1])))(x.split()) for x in los[:-1]]
        else:
            loq = [(lambda x: (units.Quantity(float(x[0]),x[1])))(x.split()) for x in los]
        try:
            _ = reduce((lambda a, b:  a.to(b.unit)), loq)
            loq = [a.to(loq[0].unit) for a in loq]
            return loq
        except UnitsException as e:
            msg = "Incompatible units in %s"%str(los) + str(e)
            raise ValueError(msg)
    
    def check_type(self, value):
        if isinstance(value, dict):
            if reduce((lambda a, b: a and b in value), [True, 'start', 'end']):
                los = [value['start'], value['end']]
                loq = self._to_units(los)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 2:
                loq = self._to_units(value)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
        return False
    
    def to_type(self, value):
        if isinstance(value, list):
            return self._to_units(value[:2])
        elif isinstance(value, dict):
            los = [value['start'], value['end']]
            return self._to_units(los)
    
    
class PropertyTypeQuantityRangeSampled(PropertyTypeQuantityRange):
    
    def check_type(self, value):
        if isinstance(value, dict):
            if reduce((lambda a, b: a and b in value), [True, 'start', 'stop', 'num']):
                los = [value['start'], value['stop']]
                loq = self._to_units(los)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 3:
                loq = self._to_units(value)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
        return False
    
    def to_type(self, value):
        if isinstance(value, list):
            _tmp = self._to_units(value[:2])
            _tmp.append(value[2])
            return _tmp
        elif isinstance(value, dict):
            los = [value['start'], value['stop']]
            _tmp = self._to_units(los)
            _tmp.append(value['num'])
            return _tmp
        
    
class PropertyTypeString(PropertyType):
    
    def check_type(self, value):
        try:
            str(value)
            if self._check_allowed_value(value):
                return True
        except ValueError:
            return False
    
    def to_type(self, value):
        return str(value)
    
class PropertyTypeStringList(PropertyTypeString):
    
    def check_type(self, value):
        try:
            str(value)
        except ValueError:
            return False
        if value in self.allowed_value:
            return True
        else:
            return False
    
    def to_type(self, value):
        return str(value)
        
        pass
    
    
class PropertyTypeList(PropertyType):
    
    def check_type(self, value):
        try:
            return isinstance(value, list)
        except ValueError:
            return False
    
    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, basestring):
            return value.split()
        else:
            return []
    
    
class PropertyTypeRange(PropertyType):
    
    def check_type(self, value):
        if isinstance(value, dict):
            if reduce((lambda a, b: a in value), [True, 'start', 'stop']):
                if abs(value['start'] - value['stop']) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 2:
                if abs(value[0] - value[1]) > 0:
                    return True
        return False
    
    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, dict):
            return [value['start'], value['stop']]
    
class PropertyTypeRangeSampled(PropertyTypeRange):
    
    def check_type(self, value):
        if isinstance(value, dict):
            if reduce((lambda a, b: a in value),\
                [True, 'start', 'stop', 'num']):
                if abs(value['start'] - value['stop']) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 3:
                if abs(value[0] - value[1]) > 0:
                    return True
        return False
        
    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, dict):
            return [value['start'], value['stop'], value['num']]
    
class PropertyTypeAbundances(PropertyType):
    
    elements = { 'neut': 0, 'h': 1, 'he': 2, 'li': 3, 'be': 4, 'b': 5, 'c': 6, 'n': 7, 'o': 8, 'f': 9, 'ne': 10, 'na': 11, 'mg': 12, 'al': 13, 'si': 14, 'p': 15, 's': 16, 'cl': 17, 'ar': 18, 'k': 19,    'ca': 20, 'sc': 21, 'ti': 22, 'v': 23, 'cr': 24, 'mn': 25, 'fe': 26, 'co': 27, 'ni': 28, 'cu': 29, 'zn': 30, 'ga': 31, 'ge': 32, 'as': 33, 'se': 34, 'br': 35, 'kr': 36, 'rb': 37, 'sr': 38, 'y': 39,  'zr': 40, 'nb': 41, 'mo': 42, 'tc': 43, 'ru': 44, 'rh': 45, 'pd': 46, 'ag': 47, 'cd': 48, 'in': 49, 'sn': 50, 'sb': 51, 'te': 52, 'i': 53, 'xe': 54, 'cs': 55, 'ba': 56, 'la': 57, 'ce': 58, 'pr': 59, 'nd': 60, 'pm': 61, 'sm': 62, 'eu': 63, 'gd': 64, 'tb': 65, 'dy': 66, 'ho': 67, 'er': 68, 'tm': 69, 'yb': 70, 'lu': 71, 'hf': 72, 'ta': 73, 'w': 74, 're': 75, 'os': 76, 'ir': 77, 'pt': 78, 'au': 79, 'hg': 80, 'tl': 81, 'pb': 82, 'bi': 83, 'po': 84, 'at': 85, 'rn': 86, 'fr': 87, 'ra': 88, 'ac': 89, 'th': 90, 'pa': 91, 'u': 92, 'np': 93, 'pu': 94, 'am': 95, 'cm': 96, 'bk': 97, 'cf': 98, 'es': 99, 'fm': 100, 'md': 101, 'no': 102, 'lr': 103, 'rf': 104, 'db': 105, 'sg': 106, 'bh': 107, 'hs': 108, 'mt': 109, 'ds':110, 'rg':111, 'cn':112 }

    def check_type(self, _value):
        value = dict((k.lower(), v) for k,v in _value.items())
        if set(value).issubset(set(self.elements)):
            return True
        else:
            return False
        
    def to_type(self, _value):
        if isinstance(_value, dict):
            value = dict((k.lower(), v) for k,v in _value.items())
            abundances = dict.fromkeys(self.elements.copy(), 0.0)
            for k in value:
                abundances[k] = value[k]
            return abundances
        else:
            raise ConfigError
        
class PropertyTypeLegacyAbundances(PropertyType):
    
    elements = { 'neut': 0, 'h': 1, 'he': 2, 'li': 3, 'be': 4, 'b': 5, 'c': 6, 'n': 7, 'o': 8, 'f': 9, 'ne': 10, 'na': 11, 'mg': 12, 'al': 13, 'si': 14, 'p': 15, 's': 16, 'cl': 17, 'ar': 18, 'k': 19,    'ca': 20, 'sc': 21, 'ti': 22, 'v': 23, 'cr': 24, 'mn': 25, 'fe': 26, 'co': 27, 'ni': 28, 'cu': 29, 'zn': 30, 'ga': 31, 'ge': 32, 'as': 33, 'se': 34, 'br': 35, 'kr': 36, 'rb': 37, 'sr': 38, 'y': 39,  'zr': 40, 'nb': 41, 'mo': 42, 'tc': 43, 'ru': 44, 'rh': 45, 'pd': 46, 'ag': 47, 'cd': 48, 'in': 49, 'sn': 50, 'sb': 51, 'te': 52, 'i': 53, 'xe': 54, 'cs': 55, 'ba': 56, 'la': 57, 'ce': 58, 'pr': 59, 'nd': 60, 'pm': 61, 'sm': 62, 'eu': 63, 'gd': 64, 'tb': 65, 'dy': 66, 'ho': 67, 'er': 68, 'tm': 69, 'yb': 70, 'lu': 71, 'hf': 72, 'ta': 73, 'w': 74, 're': 75, 'os': 76, 'ir': 77, 'pt': 78, 'au': 79, 'hg': 80, 'tl': 81, 'pb': 82, 'bi': 83, 'po': 84, 'at': 85, 'rn': 86, 'fr': 87, 'ra': 88, 'ac': 89, 'th': 90, 'pa': 91, 'u': 92, 'np': 93, 'pu': 94, 'am': 95, 'cm': 96, 'bk': 97, 'cf': 98, 'es': 99, 'fm': 100, 'md': 101, 'no': 102, 'lr': 103, 'rf': 104, 'db': 105, 'sg': 106, 'bh': 107, 'hs': 108, 'mt': 109, 'ds':110, 'rg':111, 'cn':112 }
    types = ['uniform']

    def check_type(self, _value):
        value = dict((k.lower(), v) for k,v in _value.items())
        print(value)
        if 'type' in value:
            if value['type'] in self.types:
                print('type is ok')
                tmp = value.copy()
                tmp.pop('type', None)
                if set(tmp).issubset(set(self.elements)):
                    return True
                else:
                    return False
        return False
        
    def to_type(self, _value):
        if isinstance(_value, dict):
            value = dict((k.lower(), v) for k,v in _value.items())
            abundances = dict.fromkeys(self.elements.copy(), 0.0)
            for k in value:
                abundances[k] = value[k]
            abundances['type'] = value['type']
            return abundances
        else:
            raise ConfigError
            
    
    


class DefaultParser(object):
    """Not invented here syndrome"""

    __check = {}
    __convert = {}
    __list_of_leaf_types = []
    __types = {}

    def __init__(self, default_dict):
        """Creates a new property object for the given config level
        :param default_dict: default configuration
        :return:
        """
        #create property type dict
        self.__types['arbitrary'] = PropertyType
        
        self.__types['int'] = PropertyTypeInt
        self.__register_leaf('int')
        
        self.__types['float'] = PropertyTypeFloat
        self.__register_leaf('float')
        
        self.__types['quantity'] = PropertyTypeQuantity
        self.__register_leaf('quantity')
        
        self.__types['quantity_range'] = PropertyTypeQuantityRange
        self.__register_leaf('quantity_range')
        
        self.__types['quantity_range_sampled'] = PropertyTypeQuantityRangeSampled
        self.__register_leaf('quantity_range_sampled')
        
        self.__types['string'] = PropertyTypeString
        self.__register_leaf('string')
        
        
        self.__types['range'] = PropertyTypeRange
        self.__register_leaf('range')
        
        self.__types['range_sampled'] = PropertyTypeRangeSampled
        self.__register_leaf('range_sampled')
        
        self.__types['list'] = PropertyTypeList
        self.__register_leaf('list')
        
        self.__types['container-declaration'] = PropertyTypeContainer
        self.__register_leaf('container-declaration')
        
        self.__types['container-property'] = PropertyTypeContainer
        self.__register_leaf('container-property')
        
        self.__types['abundance_set'] = PropertyTypeAbundances
        self.__register_leaf('abundance_set')
        
        self.__types['legacy-abundances'] =PropertyTypeLegacyAbundances
        self.__register_leaf('legacy-abundances')
        
        self.__mandatory = False
        self.__default_value = None

        self.__allowed_value = None
        self.__allowed_type = None
        self.__config_value = None
        self.__path = None

        self.__default_dict = default_dict

        if not 'property_type' in default_dict:
            self.__property_type = 'arbitrary'
        else:
            self.__property_type = default_dict['property_type']
            #print(self.__property_type)
            #print(self.__types.keys())
            if not self.__property_type in self.__types.keys():
                raise ValueError
        self.__type = self.__types[self.__property_type]()            

        if 'allowed_value' in default_dict:
            self.__type.allowed_value = default_dict['allowed_value']

        if 'allowed_type' in default_dict:
            self.__type.allowed_type = default_dict['allowed_type']
            
#ToDo: move all to classes
        if 'default' in default_dict:
            if default_dict['default'] != None and not default_dict['default'] in ['None','']:
                self.__type.default = default_dict['default']

        if 'mandatory' in default_dict:
            self.__type.mandatory = default_dict['mandatory']

        self.is_leaf = self.__is_leaf(self.__property_type)

    def get_default(self):
        """Returns the default value of this property, if specified.
        :return: default value
        """
        return self.__type.default

    def set_default(self, value):
        """
        Set a new default value.
        :param value: new default value
        """
        if value != None:
            if self.__type.check_type(value):
                self.__type.default = value
            else:
                raise ValueError('Default value violates property constraint. Check %s : %s' %(self.get_path_in_dict(), self.__property_type))
        else:
            self.__type.default = None

    def is_mandatory(self):
        """
        Returns True if this property is a mandatory.
        :return: mandatory
        """
        return self.__type.mandatory

    def has_default(self):
        """
        Returns True if this property has a default value
        :return: has a default value
        """
        try:
            if self.__type.default != None:
                return True
            else:
                return False
        except NameError:
            pass

    def set_path_in_dic(self, path):
        """
        Set the path to this property in the config
        :param path: path(chain of keys)
        :return:
        """
        self.__path = path

    def get_path_in_dict(self):
        """
        Returns the path of this property in the config
        :return: path
        """
        return self.__path

    def set_config_value(self, value):
        """
        Set a new value
        :param value:
        :return:
        """
        self.__config_value = value

    def get_value(self):
        """
        Returns the configuration value from the configuration.
        If the value specified in the configuration is invalid
        the default value is returned
        :return: value
        """
        if (self.__config_value is not None and
            self.__type.check_type(self.__config_value)):
            return self.__type.to_type(self.__config_value)
        else:
            if self.has_default():
                logger.warning("Value <%s> specified in the configuration violates a constraint\
                      given in the default configuration. Expected type: %s. Using the default value."%(str(self.__config_value),str(self.__property_type)))
                return self.__type.default
            else:
                raise ValueError('No default value given.')

    def is_container(self):
        """
        Returns True if this property is of type container.
        :return:
        """
        return self.__is_container()

    def get_container_dic(self):
        """
        If this property is a container it returns the corresponding
        container dictionary
        :return: container dictionary
        """
        if self.__is_container():
            return self.__container_dic

    @classmethod
    def update_container_dic(cls, container_dic, current_entry_name):
        if reduce(lambda a, b: a or b,\
                  [container_dic.has_key(i) for i in ['and', 'or']], True):
            if 'or' in container_dic:
                if current_entry_name in container_dic['or']:
                    container_dic['or'] = []
                    return container_dic
            if 'and' in container_dic:
                if current_entry_name in container_dic['and']:
                    current_entry_name['and'].remove(current_entry_name)
                    return container_dic


    def is_valid(self, value):
        if not self.__check[self.__property_type](self, value):
            return False
        if self.__allowed_value:
            if not self.__is_allowed_value(value, self.__allowed_value):
                return False
        if self.__allowed_type:
            if not self.__check_value(value, self.__lower, self.__upper):
                return False
        return True

    def __register_leaf(self, type_name):
        #print(type_name)
        if not type_name in self.__list_of_leaf_types:
            self.__list_of_leaf_types.append(type_name)

    def __is_leaf(self, type_name):
        return type_name in self.__list_of_leaf_types

    def __is_container(self):
        if self.__property_type == 'container-property':
            try:
                self.__container_dic = self.__default_dict['type']['containers']
                return True
            except KeyError:
                return False
        else:
            return False

#    __check['container-property'] = __is_container

    def __is_container_declaration(self, value):
        pass


class Container(DefaultParser):
    def __init__(self, container_default_dict, container_dict):
        """Creates a new container object.
        :param container_default_dict: Dictionary containing the default properties of the container.
        :param container_dict: Dictionary containing the configuration of the container.
        """

        #self.__register_leaf('list')
        #self.__register_leaf('int')
        #self.__register_leaf('float')
        #self.__register_leaf('quantity')
        #self.__register_leaf('string')

        self.__type = None
        self.__allowed_value = None
        self.__allowed_type = None
        self.__config_value = None
        self.__path = None

        self.__property_type = 'container-property'

        self.__default_container = {}
        self.__config_container = {}

        #check if it is a valid default container
        if not 'type' in container_default_dict:
            raise ValueError('The given default container is no valid')

        #set allowed containers
        try:
            self.__allowed_container = container_default_dict['type']['containers']
        except KeyError:
            raise ValueError('No container names specified')

        #check if the specified container in the config is allowed
        try:
            #print(container_dict)
            if not container_dict['type'] in self.__allowed_container:

                raise ValueError('Wrong container type')
            else:
                type_dict = container_dict['type']
                self.__type = container_dict['type']
        except KeyError:
            raise ValueError('No container type specified')

        #get selected container from conf
        try:
            self.__selected_container = container_dict['type']
        except KeyError:
            self.__selected_container = None
            raise ValueError('No container type specified in config')

        #look for necessary items
        entry_name = '_' + self.__selected_container
        try:
            necessary_items = container_default_dict['type'][entry_name]
        except KeyError:
            raise ValueError('Container insufficient specified')

        def parse_container_items(top_default, top_config, level_name, path):
            """Recursive parser for the container default dictionary and the container configuration dictionary.

            :param top_default: container default dictionary of the upper recursion level
            :param top_config: container configuration dictionary of the upper recursion level
            :param level_name: name(key) of the of the upper recursion level
            :param path: path in the nested container dictionary from the main level to the current level
            :return: If the current recursion level is not a leaf, the function returns a dictionary with itself for
            each branch. If the  current recursion level is a leaf the configured value and a configuration object is
            returned
            """
            tmp_conf_ob = {}
            tmp_conf_val = {}
            if isinstance(top_default, dict):
                default_property = DefaultParser(top_default)
                if default_property.is_container():
                    container_conf = get_value_by_path(top_config, path)
                    ccontainer = Container(top_default, container_conf)
                    return ccontainer.get_container_ob(), ccontainer.get_container_conf()
                elif not default_property.is_leaf:
 #                   print(top_default.items())
                    for k, v in top_default.items():
                        tmp_conf_ob[k], tmp_conf_val[k] = parse_container_items(v, top_config, k, path + [k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
                        conf_value = get_value_by_path(top_config, path)
                    except KeyError:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)

                    return default_property, default_property.get_value()

        def get_value_by_path(dict, path):
            """
            Value from a nested dictionary specified by its path.
            :param dict: nested source dictionary
            :param path: path (composed of keys) in the dictionary
            :return:
            """
            for key in path:
                dict = dict[key]
            return dict

        for item in necessary_items:
            if not item in container_dict.keys():
                raise ValueError('Entry %s is missing in container' % str(item))
            self.__default_container, self.__config_container = parse_container_items(container_default_dict[item],
                                                                                      container_dict[item], item, [])
            pass
            #go through  all items and create an conf object thereby check the conf

        self.__container_ob = self.__default_container
        self.__conf = self.__config_container

    def get_container_ob(self):
        """
        Return the container configuration object
        :return:
        """
        return self.__container_ob

    def get_container_conf(self):
        """
        Return the configuration
        :return:
        """
        self.__conf['type'] = self.__type
        return self.__conf


class Config(object):
    """
    An configuration object represents the parsed configuration.
    """


    def __init__(self, default_configuration, input_configuration):
        """Creates the configuration object.
        :param default_configuration: Default configuration dictionary
        :param input_configuration: Configuration dictionary
        """
        self.__conf_o = None
        self.__conf_v = None
        self.mandatories = {}
        self.fulfilled = {}
        self.__create_default_conf(default_configuration)
        self.__parse_config(default_configuration, input_configuration)

    @classmethod
    def from_yaml(cls, fname_config, fname_default):
        with open(fname_config) as f:
            conff = f.read()
            conf = yaml.safe_load(conff)
        with open(fname_default) as f:
            defaf = f.read()
            defa = yaml.safe_load(defaf)
        return cls(defa, conf)
    

    def __mandatory_key(self, path):
        """Return the key string for dictionary of mandatory entries
        :param path: path (composed of keys) in the dictionary
        :return: corresponding key
        """
        return ':'.join(path)
    
    

    def register_mandatory(self, name, path):
        """Register a mandatory entry
        :param name: name of the mandatory entry to be registered
        :param path: path (composed of keys) in the dictionary
        """
        self.mandatories[self.__mandatory_key(path)] = name

    def deregister_mandatory(self, name, path):
        """Register a deregistered mandatory entry
        :param name: name of the mandatory entry to be deregistered
        :param path: path (composed of keys) in the dictionary
        """
        self.fulfilled[self.__mandatory_key(path)] = name

    def is_mandatory_fulfilled(self):
        """
        Check if all mandatory entries are deregistered.
        """
        if len(set(self.mandatories.keys()) - set(self.fulfilled.keys())) <= 0:
            return True
        else:
            return False

    def __parse_config(self, default_configuration, configuration):
        """Parser for the default dictionary and the configuration dictionary.
        :param default_configuration: Default configuration dictionary
        :param configuration:  Configuration dictionary
        """

        def find_item(dict, key):
            """
            Returns the value for a specific key in a nested dictionary
            :param dict: nested dictionary
            :param key:
            """
            if key in dict: return dict[key]
            for k, v in dict.items():
                if isinstance(v, dict):
                    item = find_item(v, key)
                    if item is not None:
                        return item


        def get_property_by_path(dict, path):
            """ Returns the value for a specific path(chain of keys) in a nested dictionary
            :param dict: nested dictionary
            :param path: chain of keys as list
            """

            for key in path:
                dict = dict[key]
            return dict

        def recursive_parser(top_default, configuration, path):
            """
            Recursive parser for the default dictionary.
            :param top_default: container default dictionary of the upper recursion level
            :param configuration:  configuration dictionary
            :param path: path in the nested container dictionary from the main level to the current level
            :return: If the current recursion level is not a leaf, the function returns a dictionary with itself for
            each branch. If the  current recursion level is a leaf, the configuration value and object are
            returned
            """
            tmp_conf_ob = {}
            tmp_conf_val = {}
            if isinstance(top_default, dict):
                default_property = DefaultParser(top_default)
                if default_property.is_mandatory():
                    self.register_mandatory(self, path)
                self.deregister_mandatory(self, path)

                if default_property.is_container():
                    container_conf = get_property_by_path(configuration, path)
                    ccontainer = Container(top_default, container_conf)
                    return ccontainer.get_container_ob(), ccontainer.get_container_conf()
                elif not default_property.is_leaf:
                    no_default = self.__check_items_in_conf(get_property_by_path(configuration, path), top_default)
                    if len(no_default) > 0:
                        logger.warning('The items %s from the configuration are not specified in the default configuration'%str(no_default))
                    for k, v in top_default.items():
                        #print('>---<')
                        #print(k)
                        tmp_conf_ob[k], tmp_conf_val[k] = recursive_parser(v, configuration, path + [k])
                        #print(tmp_conf_val[k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
#                        print('get_property_by_path')
                        #print(path)
                        conf_value = get_property_by_path(configuration, path)
#                        print(conf_value)
#                        print('End:get_property_by_path')
                    except KeyError:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)
                    return default_property, default_property.get_value()


        self.__conf_o, self.__conf_v = recursive_parser(default_configuration, configuration, [])
 #       print('|\|\|\|')
 #       print(self.__conf_v)
        
    def __check_items_in_conf(self, config_dict, default_dict):
        return list(set(config_dict.keys()) - set(default_dict.keys()))
        

    def __create_default_conf(self, default_conf):
        """Returns the default configuration values as dictionary.
        :param default_conf: default configuration dictionary
        :return: default configuration values
        """

        def recursive_default_parser(top_default, path):
            """Recursive parser for the default dictionary.
            :param top_default: container default dictionary of the upper recursion level
            :param path: path in the nested container dictionary from the main level to the current level
            :return: If the current recursion level is not a leaf, the function returns a dictionary with itself for
            each branch. If the  current recursion level is a leaf, the default configuration value is
            returned
            """
            tmp_default = {}
            if isinstance(top_default, dict):
                default_property = DefaultParser(top_default)
                if not default_property.is_container():
                    if not default_property.is_leaf:
                        for k, v in top_default.items():
                            tmp_default[k] = recursive_default_parser(v, path + [k])
                        return tmp_default
                    else:
                        default_property.set_path_in_dic(path)
                        if default_property.has_default():
                            return default_property.get_default()
                        else:
                            return None

        self.__default_config = recursive_default_parser(default_conf, [])

    def get_config(self):
        """Returns the parsed configuration as dictionary.
        :return: configuration values as dictionary
        """
        print(self.__conf_v)
        return self.__conf_v

    def get_default_config(self):
        """Returns the default configuration values as dictionary
        :return: default configuration values as dictionary
        """
        return self.__default_config

    def get_config_object(self):
        """Returns the default configuration objects as dictionary
        :return: default configuration objects as dictionary
        """
        return self.__conf_o

