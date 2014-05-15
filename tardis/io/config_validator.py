# coding=utf-8

import re
import logging
import pprint
import ast

import numpy as np
from astropy import units

try:
    from astropy.units.core import UnitsException
except ImportError:
    from astropy.units.core import UnitsError as UnitsException

from astropy import constants

import yaml

from tardis.atomic import atomic_symbols_data


logger = logging.getLogger(__name__)


class Error(Exception):
    """Base class for exceptions in the config parser."""
    pass


class ConfigTypeError(Error, ValueError):
    """
    Exception raised if the type of the configured value mismatches the type
    specified in the default configuration.
    """

    def __init__(self, value, expected_type, _help):
        self.value = value
        self.expected_type = expected_type
        self.help = _help

    def __str__(self):
        return "Expected type %s but found %s.\nHelp:%s " % \
               (repr(self.expected_type), repr(type(self.value)), help)


class ConfigError(Error):
    """
    Exception raised if something is wrong in the configuration.
    """

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return "Error in the configuration at %s " % ("->".join(self.path))


class ConfigValueError(ConfigError, ValueError):
    """
    Exception raised if the given value does not match the allowed constraints.
    """

    default_msg = "Given value (%s) not allowed in constraint (%s). [%s]"

    def __init__(self, config_value, allowed_constraint, path, msg=None):
        self.config_value = config_value
        self.allowed_constraint = allowed_constraint
        self.path = path
        if msg is None:
            self.msg = self.default_msg
        else:
            self.msg = msg

    def __str__(self):
        return self.msg % (str(self.config_value), str(self.allowed_constraint), self.path)


class DefaultConfigError(ConfigError):
    """
    Exception raised if something is wrong in the default configuration.
    """

    def __str__(self):
        return "Error in the default configuration at %s " % \
               ("->".join(self.path))


class PropertyType(object):
    """
    Base class for all property types containing all the basic methods.
    """

    def __init__(self):
        self._default = None
        self._allowed_value = None
        self._allowed_type = None
        self._help = None
        self._mandatory = False
        self._lower = None
        self._upper = None
        pass

    @property
    def default(self):
        """
        Geter for the default config value.
        Returns
        -------
        default
                default config value.
        """
        return self._default

    @default.setter
    def default(self, value):
        """
        Sets the default value if the type is ok.
        Parameters
        ----------
        value
                default value
        """
        self._default = self.to_type(value)

    @property
    def allowed_value(self):
        """
        Returns the allowed value
        Returns
        -------
        allowed_value
                        allowed value
        """
        return self._allowed_value

    @allowed_value.setter
    def allowed_value(self, value):
        """
        Sets the allowed values
        Parameters
        ----------
        value
                allowed values
        """
        if isinstance(value, basestring):
            self._allowed_value = set(self.__list_dtype(value.split()))
        elif isinstance(value, list) or isinstance(value, set):
            self._allowed_value = set(self.__list_dtype(value))
        elif isinstance(value, float) or isinstance(value, int):
            self._allowed_type = {value}
        else:
            raise ValueError("Can not set allowed value.")

    @property
    def allowed_type(self):
        """
        Returns the allowed type
        Returns
        -------
        allowed_type
                        allowed type

        """
        return self._allowed_value

    @allowed_type.setter
    def allowed_type(self, value):
        self._allowed_type = value
        if '_parse_allowed_type' in (set(dir(self.__class__)) - set(dir(PropertyType))) and value is not None:
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

    @staticmethod
    def __list_dtype(mixed_list):
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
        if self._allowed_value is not None:
            atype = type(iter(self.allowed_value).next())
            value = atype(_value)
            if value in self.allowed_value:
                return True
            else:
                return False
        else:
            return True

    def __repr__(self):
        if hasattr(self, "_allowed_type"):
            return str("Type %s; Allowed type: %s" % (self.__class__.__name__, self._allowed_type))
        else:
            return str("Type %s; " % self.__class__.__name__)


class PropertyTypeContainer(PropertyType):
    def check_type(self):
        pass


class PropertyTypeBool(PropertyType):
    def check_type(self, value):
        try:
            foo = bool(value)
            return True
        except ValueError:
            return False

    def to_type(self, value):
        return bool(value)


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

    @staticmethod
    def _parse_allowed_type(allowed_type):
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

    @staticmethod
    def __check_type(value, lower_lim, upper_lim):
        upper, lower = True, True
        if upper_lim is not None:
            upper = value < upper_lim
        if lower_lim is not None:
            lower = value > lower_lim
        return upper and lower

    def _check_allowed_type(self, value):
        if self._allowed_type is not None:
            if self.__check_type(value, self._lower, self._upper):
                return True
            else:
                return False
        else:
            return True

    def _is_valid(self, value):
        if not self.check_type(value):
            return False
        if self.allowed_value is not None:
            return False
        if not self.__check_type(value, self._lower, self._upper):
            return False
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
            quantity_split = value.strip().split()
            quantity_value = quantity_split[0]
            quantity_unit = ' '.join(quantity_split[1:])
            try:
                float(quantity_value)
                units.Unit(quantity_unit)
            except ValueError:
                return False

            if self._default is not None:
                #d_quantity_split = self._default.strip().split()
                self._default.to(quantity_unit)
            float(quantity_value)
            units.Unit(quantity_unit)
            return True
        except (ValueError, AttributeError):
            return False

    def to_type(self, value):
        quantity_split = value.strip().split()
        quantity_value = quantity_split[0]
        quantity_unit = ' '.join(quantity_split[1:])
        if quantity_unit.strip() == 'log_lsun':
            quantity_value = 10 ** (float(quantity_value) +
                                    np.log10(constants.L_sun.cgs.value))
            quantity_unit = 'erg/s'

        return float(quantity_value) * units.Unit(quantity_unit)


class PropertyTypeQuantityRange(PropertyTypeQuantity):
    @staticmethod
    def _to_units(los):
        if len(los) > 2:
            loq = [(lambda x: (units.Quantity(float(x[0]), x[1])))(x.split())
                   for x in los[:-1]]
        else:
            loq = [(lambda x: (units.Quantity(float(x[0]), x[1])))(x.split())
                   for x in los]
        try:
            _ = reduce((lambda a, b: a.to(b.unit)), loq)
            loq = [a.to(loq[0].unit) for a in loq]
            return loq
        except UnitsException as e:
            msg = "Incompatible units in %s" % str(los) + str(e)
            raise ValueError(msg)

    def check_type(self, value):
        if isinstance(value, dict):
            if (reduce((lambda a, b: a and b in value.keys()), [True, 'start', 'end'])) \
                    or (reduce((lambda a, b: a and b in value.keys()), [True, 'start', 'stop'])):  # for legacy support
                los = [value['start'], value['end']]
                loq = self._to_units(los)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 2:
                loq = self._to_units(value)
                if abs(loq[0].value - loq[1].value) > 0:
                    return True
                else:
                    return False
        elif isinstance(value, basestring):
            try:
                clist = ast.literal_eval(value)
                if len(clist) == 2:
                    loq = self._to_units(value)
                    if abs(loq[0].value - loq[1].value) > 0:
                        return True
                    else:
                        return False
                else:
                    return False
            except SyntaxError:
                clist = value.split()
                if len(clist) == 2:
                    loq = self._to_units(value)
                    if abs(loq[0].value - loq[1].value) > 0:
                        return True
                    else:
                        return False
                else:
                    return False
            except ValueError:
                return False
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
        if isinstance(value, list):
            return True
        elif isinstance(value, basestring):
            try:
                ast.literal_eval(value)
                return True
            except SyntaxError:
                try:
                    value.split()
                    return True
                except AttributeError:
                    return False
        return False

    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, basestring):
            try:
                return ast.literal_eval(value)
            except SyntaxError:
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
        elif isinstance(value, basestring):
            try:
                clist = ast.literal_eval(value)
                if abs(clist[0] - clist[1]) > 0:
                    return True
            except SyntaxError:
                clist = value.split()
                if abs(clist[0] - clist[1]) > 0:
                    return True
        return False

    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, dict):
            return [value['start'], value['stop']]
        elif isinstance(value, basestring):
            try:
                return ast.literal_eval(value)
            except SyntaxError:
                return value.split()


class PropertyTypeRangeSampled(PropertyTypeRange):
    def check_type(self, value):
        if isinstance(value, dict):
            if reduce((lambda a, b: a in value),
                      [True, 'start', 'stop', 'num']):
                if abs(value['start'] - value['stop']) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 3:
                if abs(value[0] - value[1]) > 0:
                    return True
        elif isinstance(value, basestring):
            try:
                clist = ast.literal_eval(value)
                if abs(clist[0] - clist[1]) > 0:
                    return True
            except SyntaxError:
                clist = value.split()
                if abs(clist[0] - clist[1]) > 0:
                    return True
        return False

    def to_type(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, dict):
            return [value['start'], value['stop'], value['num']]
        elif isinstance(value, basestring):
            try:
                return ast.literal_eval(value)
            except SyntaxError:
                return value.split()


class PropertyTypeAbundances(PropertyType):
    elements = dict([(x, y.lower()) for (x, y) in atomic_symbols_data])

    def check_type(self, _value):
        if isinstance(_value, dict):
            value = dict((k.lower(), v) for k, v in _value.items())
            if set(value).issubset(set(self.elements.values())):
                return True
            else:
                return False
        else:
            return False

    def to_type(self, _value):
        if isinstance(_value, dict):
            value = dict((k.lower(), v) for k, v in _value.items())
            abundances = dict.fromkeys(self.elements.copy(), 0.0)
            for k in value:
                abundances[k] = value[k]
                abundances = {k: v for k, v in abundances.items() if not (v == 0.)}
            return abundances
        else:
            raise ConfigError


class PropertyTypeLegacyAbundances(PropertyType):
    elements = dict([(x, y.lower()) for (x, y) in atomic_symbols_data])
    types = ['uniform']

    def check_type(self, _value):
        value = dict((k.lower(), v) for k, v in _value.items())
        if 'type' in value:
            if value['type'] in self.types:
                tmp = value.copy()
                tmp.pop('type', None)
                if set(tmp).issubset(set(self.elements.values())):
                    return True
                else:
                    return False
        return False

    def to_type(self, _value):
        if isinstance(_value, dict):
            value = dict((k.lower(), v) for k, v in _value.items())
            abundances = dict.fromkeys(self.elements.copy(), 0.0)
            for k in value:
                abundances[k] = value[k]
            abundances['type'] = value['type']
            return abundances
        else:
            raise ConfigError


class DefaultParser(object):
    """
    Creates a new property object for the given config level

    Parameters
    ----------
    default_dict: dict
        default configuration


    """

    __check = {}
    __convert = {}
    __list_of_leaf_types = []
    __types = {}

    def __init__(self, default_dict, item_path=None):

        self.__item_path = item_path
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

        self.__types['legacy-abundances'] = PropertyTypeLegacyAbundances
        self.__register_leaf('legacy-abundances')

        self.__types['bool'] = PropertyTypeBool
        self.__register_leaf('bool')

        self.__mandatory = False

        self.__default_value = None

        self.__allowed_value = None
        self.__allowed_type = None
        self.__config_value = None
        self.__path = None

        self.__container_dic = None

        self.__default_dict = default_dict

        if not 'property_type' in default_dict:
            self.__property_type = 'arbitrary'
        else:
            self.__property_type = default_dict['property_type']
            if not self.__property_type in self.__types.keys():
                raise ValueError
        self.__type = self.__types[self.__property_type]()

        if 'allowed_value' in default_dict:
            self.__type.allowed_value = default_dict['allowed_value']

        if 'allowed_type' in default_dict:
            self.__type.allowed_type = default_dict['allowed_type']

        if 'default' in default_dict:
            if default_dict['default'] is not None and not default_dict['default'] in ['None', '']:
                self.__type.default = default_dict['default']

        if 'mandatory' in default_dict:
            self.__type.mandatory = default_dict['mandatory']

        self.is_leaf = self.__is_leaf(self.__property_type)

    def get_default(self):
        """Returns the default value of this property, if specified.

        Returns
        -------
        default value
        """
        return self.__type.default

    def set_default(self, value):
        """Set a new default value.
        Parameters
        ----------
        value:
                new default value
        """
        if value is not None:
            if self.__type.check_type(value):
                self.__type.default = value
            else:
                raise ConfigValueError(value, self.__type.allowed_value, self.get_path_in_dict(),
                                       msg='Default value (%s) violates property constraint (%s). [%s]')
        else:
            self.__type.default = None

    @property
    def is_mandatory(self):
        """Returns True if this property is a mandatory.
        Returns
        -------
        bool
                True if this property is a mandatory.
        """
        return self.__type.mandatory

    @property
    def has_default(self):
        """Returns True if this property has a default value
        Returns
        -------
        bool
                True if a default value was given.
        """
        try:
            if self.__type.default is not None:
                return True
            else:
                return False
        except NameError:
            pass

    def set_path_in_dic(self, path):
        """Set the path to this property in the config
        Parameters
        ----------
        path:   list  of str
                Path in config dictionary.
        """
        self.__path = path

    def get_path_in_dict(self):
        """Returns the path of this property in the config
        Returns
        -------
        path:   list of str
                Path in config dictionary.
        """
        return self.__path

    def set_config_value(self, value):
        """Set a new config value.
        Parameters
        ----------
        value
                New config value.
        """
        self.__config_value = value

    def get_value(self):
        """
        Returns the configuration value from the configuration.
        If the value specified in the configuration is invalid
        the default value is returned
        Returns
        -------
        value
                Config value.
        """
        if self.__config_value is not None:
            if self.__type.check_type(self.__config_value):
                return self.__type.to_type(self.__config_value)
            else:
                raise ConfigValueError(self.__config_value, self.__type.allowed_value, self.get_path_in_dict())
        else:
            if self.has_default:
                logger.debug("Value <%s> specified in the configuration violates a constraint "
                             "given in the default configuration. Expected type: %s. Using the default value." % (
                                 str(self.__config_value), str(self.__property_type)))
                return self.__type.default
            else:
                if self.is_mandatory:
                    raise ValueError('Value is mandatory, but no value was given in default configuration. [%s]' % str(
                        self.get_path_in_dict()))
                else:
                    logger.debug("Value is not mandatory and is not specified in the configuration. [%s]" % (
                        str(self.get_path_in_dict())))
                    return None

    def is_container(self):
        """
        Returns True if this property is of type container.
        """
        return self.__is_container()

    def get_container_dic(self):
        """
        If this property is a container it returns the corresponding
        container dictionary
        Returns
        -------
        container dictionary:   dict
                                Container dictionary
        """
        if self.__is_container():
            return self.__container_dic

    @classmethod
    def update_container_dic(cls, container_dic, current_entry_name):
        if reduce(lambda a, b: a or b,
                  [(i in container_dic) for i in ['and', 'or']], True):
            if 'or' in container_dic:
                if current_entry_name in container_dic['or']:
                    container_dic['or'] = []
                    return container_dic
            if 'and' in container_dic:
                if current_entry_name in container_dic['and']:
                    current_entry_name['and'].remove(current_entry_name)
                    return container_dic

    def __register_leaf(self, type_name):
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
    def __init__(self, container_default_dict, container_dict, container_path=None):
        """Creates a new container object.
        Parameters
        ----------
        container_default_dict: dict
                                Dictionary containing the default properties of the container.
        container_dict:         dict
                                Dictionary containing the configuration of the container.
        """

        self.__container_path = container_path
        self.__type = None
        self.__allowed_value = None
        self.__allowed_type = None
        self.__config_value = None
        self.__path = None
        self.__paper_abundances = False
        self.__has_additional_items = False

        self.__property_type = 'container-property'

        self.__default_container = {}
        self.__config_container = {}

        #check if it is a valid default container
        if not 'type' in container_default_dict:
            raise ValueError('No type specified in the default configuration. [%s]' % self.__container_path)

        #set allowed containers
        try:
            self.__allowed_container = container_default_dict['type']['containers']
        except KeyError:
            raise ValueError('No container names specified in the default configuration. [%s]' % self.__container_path)

        #check if the specified container is given in the config
        if container_dict is None:
            logger.debug('%s specified in the default configuration but not present in the configuration. [%s]'
                         % (container_path[-1], container_path))
            self.__container_ob = None
            self.__conf = None
            return
        #check if the specified container in the config is allowed
        try:
            if not container_dict['type'] in self.__allowed_container:

                raise ValueError(
                    'Wrong container type in the configuration! The container is %s but only %s'
                    ' are allowed containers. [%s]' % (container_dict['type'], self.__allowed_container,
                                                       self.__container_path))
            else:
                type_dict = container_dict['type']
                self.__type = container_dict['type']
        except KeyError:
            raise ValueError('No type specified in the configuration. [%s]' % self.__container_path)

        #get selected container from conf
        try:
            self.__selected_container = container_dict['type']
        except KeyError:
            self.__selected_container = None
            raise ValueError('No container type specified in config')

        ####This is for the uniform abundances section in the paper.
        if self.__type == 'uniform' and self.__container_path[-1] == 'abundances':
            logger.debug('Using the legacy support for the uniform abundances section in the paper.')
            self.__paper_abundances = True
            cabundances_section = PropertyTypeAbundances()
            tmp_container_dict = dict(container_dict)
            tmp_container_dict.pop('type', None)
            cabundances_section.check_type(tmp_container_dict)
            tmp = cabundances_section.to_type(tmp_container_dict)
            self.__default_container, self.__config_container = tmp, tmp
        ####
        else:
            #look for necessary items
            entry_name = '_' + self.__selected_container
            try:
                necessary_items = container_default_dict['type'][entry_name]
            except KeyError:
                raise ValueError('Container insufficient specified in the default configuration. The declaration of'
                                 ' necessary items is missing. Add %s: ["your items"] to %s' % (
                    necessary_items, self.__container_path))

                #look for additional items
            entry_name = '+' + self.__selected_container
            self.__has_additional_items = False
            try:
                additional_items = container_default_dict['type'][entry_name]
                self.__has_additional_items = True
            except KeyError:
                self.__has_additional_items = False
                pass

        if not self.__paper_abundances:
            for nitem in necessary_items:
                if not nitem in container_dict:
                    raise ValueError('Entry %s is missing in configuration. [%s]' % (str(nitem), self.__container_path))
                else:
                    self.__default_container[nitem], self.__config_container[nitem] = self.parse_container_items(
                        container_default_dict[nitem],
                        container_dict[nitem], nitem, self.__container_path + [nitem])
            if self.__has_additional_items:
                for aitem in additional_items:
                    try:
                        if aitem in container_dict:
                            self.__default_container[aitem], self.__config_container[aitem] = \
                                self.parse_container_items(container_default_dict[aitem],
                                                           container_dict[aitem], aitem,
                                                       self.__container_path + [aitem])
                        else:
                            self.__default_container[aitem], self.__config_container[aitem] = \
                                self.parse_container_items(container_default_dict[aitem],
                                                           None, aitem,
                                                           self.__container_path + [aitem])

                    except KeyError:
                        pass

                            #go through  all items and create an conf object thereby check the conf
        self.__container_ob = self.__default_container
        if isinstance(self.__config_container, dict):
            self.__conf = self.__config_container
        else:
            self.__conf = {"No Name": self.__config_container}

    def parse_container_items(self, top_default, top_config, item, full_path):
        """Recursive parser for the container default dictionary and the container configuration
        dictionary.

        Parameters
        ----------
        top_default:    dict
                        container default dictionary of the upper recursion level
        top_config:     dict
                        container configuration dictionary of the upper recursion level
        level_name:     str
                        name(key) of the of the upper recursion level
        path:           list of str
                        path in the nested container dictionary from the main level to the current level
        Returns
        -------
                        If the current recursion level is not a leaf, the function returns a dictionary with itself for
        each branch. If the  current recursion level is a leaf the configured value and a configuration object is
        returned
        """

        def reduce_list(a, b):
            """
            removes items from list a which are in b. (o.B.d.A trivial)
            Parameters
            ----------
            a:  list
                minuend
            b:  list
                subtrahend
            Returns
            -------
            a:  list
                difference
            """
            for k in b:
                a.remove(k)
            return a

        def get_value_by_path(_dict, path):
            """
            Value from a nested dictionary specified by its path.
            Parameters
            ----------
            dict:   dict
                    nested source dictionary
            path:   list of str
                    path (composed of keys) in the dictionary
            Returns
            -------
            dict:   str
                    value corresponding to the given path
            """
            if _dict is None:
                return None
            else:
                for key in path:
                    _dict = _dict[key]
                return _dict

        path = reduce_list(list(full_path), self.__container_path + [item])
        tmp_conf_ob = {}
        tmp_conf_val = {}
        if isinstance(top_default, dict):
            default_property = DefaultParser(top_default)
            if default_property.is_container():
                container_conf = get_value_by_path(top_config, path)
                ccontainer = Container(top_default, container_conf, container_path=full_path)
                return ccontainer.get_container_ob(), ccontainer.get_container_conf()
            elif not default_property.is_leaf:
                for k, v in top_default.items():
                    if top_config is None:
                        tmp_conf_ob[k], tmp_conf_val[k] = None, None
                    else:
                        self.__container_path = list(full_path)
                        tmp_conf_ob[k], tmp_conf_val[k] = self.parse_container_items(v, top_config[k], k,
                                                                                     full_path + [k])

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

    def get_container_ob(self):
        """
        Return the container configuration object
        Returns
        -------
        self.__container_ob:    DefaultParser
                                container configuration object
        """
        return self.__container_ob

    def get_container_conf(self):
        """
        Return the configuration
        Returns
        -------
        self.__container_ob:    dict
                                container configuration
        """
        if self.__conf is not None:
            self.__conf['type'] = self.__type
            return self.__conf
        else:
            return self.__conf


class Config(object):
    """
    An configuration object represents the parsed configuration.

    Creates the configuration object.
    Parameters
    ----------
    configuration_definition:  dict
                            Default configuration dictionary
    input_configuration:    dict
                            Configuration dictionary
    """

    def __init__(self, configuration_definition, input_configuration):
        self.__conf_o = None
        self.__conf_v = None
        self.mandatories = {}
        self.fulfilled = {}
        self.__create_default_conf(configuration_definition)
        self.__parse_config(configuration_definition, input_configuration)
        self.__help_string = ''
        self.__default_config = None

    @classmethod
    def from_yaml(cls, fname_config, fname_default):
        with open(fname_config) as f:
            conff = f.read()
            conf = yaml.safe_load(conff)
        with open(fname_default) as f:
            defaf = f.read()
            defa = yaml.safe_load(defaf)
        return cls(defa, conf)

    @staticmethod
    def __mandatory_key(path):
        """Return the key string for dictionary of mandatory entries
        Parameters
        ----------
        path:           list of str
                        path (composed of keys) in the dictionary
        Returns
        -------
        mandatory_key:  str
                        corresponding key
        """
        return ':'.join(path)

    def register_mandatory(self, name, path):
        """Register a mandatory entry
        Parameters
        ----------
        name:   str
                name of the mandatory entry to be registered
        path:   list of str
                path (composed of keys) in the dictionary
        """
        self.mandatories[self.__mandatory_key(path)] = name

    def deregister_mandatory(self, name, path):
        """Register a deregistered mandatory entry
        Parameters
        ----------
        name:   str
                name of the mandatory entry to be deregistered
        path:   list of str
                path (composed of keys) in the dictionary
        """
        self.fulfilled[self.__mandatory_key(path)] = name

    def is_mandatory_fulfilled(self):
        """
        Check if all mandatory entries are deregistered.
        Returns
        -------
        mandatory:  bool
                    True if all mandatory entries are deregistered, otherwise False
        """
        if len(set(self.mandatories.keys()) - set(self.fulfilled.keys())) <= 0:
            return True
        else:
            return False

    def __parse_config(self, default_configuration, configuration):
        """Parser for the default dictionary and the configuration dictionary.
        Parameters
        ------
        default_configuration:  dict
                                Default configuration dictionary
        configuration:          dict
                                Configuration dictionary
        """

        def find_item(_dict, key):
            """
            Returns the value for a specific key in a nested dictionary
            note:: Deprecated, use get_property_by_path()
            Parameters
            ----------
            _dict:   dict
                    nested dictionary
            key:    str
                    key in the nested dictionary
            Returns
            -------
            item:   object
                    value corresponding to the specific key.
            """
            if key in _dict:
                return _dict[key]
            for k, v in _dict.items():
                if isinstance(v, _dict):
                    item = find_item(v, key)
                    if item is not None:
                        return item

        def get_property_by_path(d, path):
            """ Returns the value for a specific path(chain of keys) in a nested dictionary
            Parameters
            ----------
            dict:   dict
                    nested dictionary
            path:   list of str
                    chain of keys as list
            Returns
            -------
            item:   object
                    value in the nested dictionary at the specific path.
            """
            if len(path) <= 0:
                return d
            else:
                try:
                    v = d
                    for k in path:
                        v = v[k]
                    return v
                except KeyError:
                    return None

        def recursive_parser(top_default, configuration, path):
            """
            Recursive parser for the default dictionary.
            Parameters
            ----------
            top_default:    dict
                            container default dictionary of the upper recursion level
            configuration:  dict
                            configuration dictionary
            path:           list of str
                            path in the nested container dictionary from the main level to the current level
            Returns
            -------
            If the current recursion level is not a leaf, the function returns a dictionary with itself for
            each branch. If the  current recursion level is a leaf, the configuration value and object are
            returned
            """
            tmp_conf_ob = {}
            tmp_conf_val = {}
            if isinstance(top_default, dict):
                default_property = DefaultParser(top_default, item_path=path)
                if default_property.is_mandatory:
                    self.register_mandatory(self, path)
                self.deregister_mandatory(self, path)

                if default_property.is_container():
                    container_conf = get_property_by_path(configuration, path)
                    #try:
                    ccontainer = Container(top_default, container_conf, container_path=path)
                    return ccontainer.get_container_ob(), ccontainer.get_container_conf()
                #ToDo: remove general except!!!
                #except:
                #    logger.warning('Container specified in default_configuration, but not used in the current\
                #configuration file. [%s]' % str(path))
                #    return None, None

                elif not default_property.is_leaf:
                    no_default = self.__check_items_in_conf(get_property_by_path(configuration, path), top_default)
                    if len(no_default) > 0:
                        logger.debug('The items %s from the configuration are not specified in the default '
                                     'configuration' % str(no_default))
                    for k, v in top_default.items():
                        tmp_conf_ob[k], tmp_conf_val[k] = recursive_parser(v, configuration, path + [k])
                    return tmp_conf_ob, tmp_conf_val

                else:
                    default_property.set_path_in_dic(path)
                    try:
                        conf_value = get_property_by_path(configuration, path)
                    except KeyError:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)
                    return default_property, default_property.get_value()

        self.__conf_o, self.__conf_v = recursive_parser(default_configuration, configuration, [])

    @staticmethod
    def __check_items_in_conf(config_dict, default_dict):
        if isinstance(config_dict, dict) and len(config_dict) > 0:
            return list(set(config_dict.keys()) - set(default_dict.keys()))
        else:
            return list(default_dict.keys())

    def __create_default_conf(self, default_conf):
        """Returns the default configuration values as dictionary.
        Parameters
        ----------
        default_conf:   dict
                        default configuration dictionary
        Returns
        -------
        default configuration values
        """

        def recursive_default_parser(top_default, path):
            """Recursive parser for the default dictionary.
            Parameters
            ----------
            top_default:    dict
                            container default dictionary of the upper recursion level
            path:           list of str
                            path in the nested container dictionary from the main level to the current level
            Returns
            -------
            If the current recursion level is not a leaf, the function returns a dictionary with itself for
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
                        if default_property.has_default:
                            return default_property.get_default()
                        else:
                            return None

        self.__default_config = recursive_default_parser(default_conf, [])

    def get_config(self):
        """Returns the parsed configuration as dictionary.
        Returns
        -------
        configuration:  dict
                        configuration values as dictionary
        """
        return self.__conf_v

    def get_default_config(self):
        """Returns the default configuration values as dictionary
        Returns
        -------
        default_configuration:  dict
                                default configuration values as dictionary
        """
        return self.__default_config

    def get_config_object(self):
        """Returns the default configuration objects as dictionary
        Returns
        -------
        default_configuration:  objects
                                default configuration objects as dictionary
        """
        return self.__conf_o

    def get_help(self):
        return pprint.pformat(self.get_default_config())

    def __repr__(self):
        return str(pprint.pformat(self.get_config()))
        

