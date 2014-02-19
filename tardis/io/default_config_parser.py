# coding=utf-8

import re
from astropy import units



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


class DefaultParser(object):
    """Not invented here syndrome"""

    __check = {}
    __convert = {}
    __list_of_leaf_types = []

    def __init__(self, default_dict):
        """Creates a new property object for the given config level
        :param default_dict: default configuration
        :return:
        """
        self.__register_leaf('list')
        self.__register_leaf('int')
        self.__register_leaf('float')
        self.__register_leaf('quantity')
        self.__register_leaf('string')
        self.__register_leaf('container-declaration')
        self.__register_leaf('range')
        self.__register_leaf('range_sampled')
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
            if not self.__property_type in self.__check:
                raise ValueError

        if 'allowed_value' in default_dict:
            self.__allowed_value = self.__convert_av_in_pt(
                default_dict['allowed_value'], self.__property_type)

        if 'allowed_type' in default_dict:
            self.__allowed_type = default_dict['allowed_type']
            self.__lower, self.__upper = self.__parse_allowed_type(
                self.__allowed_type)

        if 'default' in default_dict:
            self.set_default(default_dict['default'])

        if 'mandatory' in default_dict:
            self.__mandatory = default_dict['mandatory']

        self.is_leaf = self.__is_leaf(self.__property_type)

    def get_default(self):
        """Returns the default value of this property, if specified.
        :return: default value
        """
        return self.__default_value

    def set_default(self, value):
        """
        Set a new default value.
        :param value: new default value
        """
        if value != None:
            if self.is_valid(value):
                self.__default_value = value
            else:
                raise ValueError('Default value violates property constraint.')
        else:
            self.__default_value = None

    def is_mandatory(self):
        """
        Returns True if this property is a mandatory.
        :return: mandatory
        """
        return self.__mandatory

    def has_default(self):
        """
        Returns True if this property has a default value
        :return: has a default value
        """
        try:
            if self.__default_value:
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
            self.is_valid(self.__config_value)):
            return self.__config_value
        else:
            if self.has_default():
                return self.__default_value
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
        print(type_name)
        if not type_name in self.__list_of_leaf_types:
            self.__list_of_leaf_types.append(type_name)

    def __is_leaf(self, type_name):
        return type_name in self.__list_of_leaf_types

    def __is_container(self):
        if self.__property_type == 'container-property':
            try:
                self.__container_dic = self.__default_dict['type']['containers']
                return True
            except:
                return False
        else:
            return False

    __check['container-property'] = __is_container

    def __is_container_declaration(self, value):
        pass


    def __is_type_arbitrary(self, value):
        self.is_leaf = False
        return True

    __check['arbitrary'] = __is_type_arbitrary

    def __is_type_list(self, value):
        self.__register_leaf('list')
        try:
            return isinstance(value, list)
        except ValueError:
            return False

    __check['list'] = __is_type_list

    def __is_type_int(self, value):
        self.__register_leaf('int')
        try:
            int(value)
            if float.is_integer(float(value)):
                return True
            else:
                return False
        except ValueError:
            return False

    __check['int'] = __is_type_int

    def __is_type_float(self, value):
        self.__register_leaf('float')
        try:
            float(value)
            return True
        except ValueError:
            return False

    __check['float'] = __is_type_float

    def __is_type_quantity(self, value):
        self.__register_leaf('quantity')
        try:
            quantity_value, quantity_unit = value.strip().split()
            float(quantity_value)
            units.Unit(quantity_unit)
            return True
        except ValueError:
            return False

    __check['quantity'] = __is_type_quantity

    def __is_type_string(self, value):
        self.__register_leaf('string')
        try:
            str(value)
            return True
        except ValueError:
            return False

    __check['string'] = __is_type_string


    def __is_type_range(self, value):
        print('----')
        print(value)
        if isinstance(value, dict):
            if reduce((lambda a, b: a in value), [True, 'start', 'end']):
                if abs(value['start'] - value['end']) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 2:
                if abs(value[0] - value[1]) > 0:
                    return True
        return False

    __check['range'] = __is_type_range


    def __is_type_range_sampled(self, value):
        print('----')
        print(value)
        if isinstance(value, dict):
            if reduce((lambda a, b: a in value),\
                [True, 'start', 'end', 'sample']):
                if abs(value['start'] - value['end']) > 0:
                    return True
        elif isinstance(value, list):
            if len(value) == 3:
                if abs(value[0] - value[1]) > 0:
                    return True
        return False

    __check['range_sampled'] = __is_type_range_sampled

    def __to_range(self, value):
        return [value['start'], value['end']]

    __convert['range'] = __to_range

    def __to_range_sampled(self, value):
        return [value['start'], value['end'], value['sample']]


    def __to_quantity(self, value):
        quantity_value, quantity_unit = value.strip().split()
        float(quantity_value)
        units.Unit(quantity_unit)
        return (quantity_value, quantity_unit)

    __convert['quantity'] = __to_quantity

    def __to_int(self, value):
        return int(value)

    __convert['int'] = __to_int

    def __to_float(self, value):
        return float(value)

    __convert['float'] = __to_float
    
    def __to_string(self, value):
        return str(value)

    __convert['string'] = __to_string

    def __to_list(self, value):
        if isinstance(value, list):
            return value
        elif isinstance(value, basestring):
            return value.split()
        else:
            return []

    __convert['list'] = __to_list

    def __convert_av_in_pt(self, allowed_value, property_type):
        """
        Converts the allowed values to the property type.
        """
        if not len([]) == 0:
            return [self.__convert[property_type](a) for a in allowed_value]
        else:
            return []

    def __is_allowed_value(self, value, allowed_value):
        if value in allowed_value:
            return True
        else:
            return False

    def __parse_allowed_type(self, allowed_type):
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

    def __check_value(self, value, lower_lim, upper_lim):

        upper, lower = True, True
        if upper_lim != None:
            upper = value < upper_lim
        if lower_lim != None:
            lower = value > lower_lim
        return upper and lower



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
        except:
            raise ValueError('No container names specified')

        #check if the specified container in the config is allowed
        try:
            if not container_dict['type'] in self.__allowed_container:

                raise ValueError('Wrong container type')
            else:
                type_dict = container_dict['type']
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
                    print(top_default.items())
                    for k, v in top_default.items():
                        tmp_conf_ob[k], tmp_conf_val[k] = parse_container_items(v, top_config, k, path + [k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
                        conf_value = get_value_by_path(top_config, path)
                    except:
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
        return self.__conf


class Config:
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
                    for k, v in top_default.items():
                        tmp_conf_ob[k], tmp_conf_val[k] = recursive_parser(v, configuration, path + [k])
                        print('>---<')
                        print(k)
                        print(tmp_conf_val[k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
                        print('get_property_by_path')
                        print(path)
                        conf_value = get_property_by_path(configuration, path)
                        print(conf_value)
                        print('End:get_property_by_path')
                    except:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)
                    return default_property, default_property.get_value()


        self.__conf_o, self.__conf_v = recursive_parser(default_configuration, configuration, [])
        print('|\|\|\|')
        print(self.__conf_v)

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


"""
class PropertyBase:
    children = {}

    def __init__(self, default_dic, config_dic, default_section, is_bottom_of_config_dict=False):
        print("")
        print("")
        self.__default_property = DefaultParser(default_dic)
        if self.__default_property.is_leaf:
            pass
            #Check conf or get help,default
        else:
            tmp = {}
            for child in default_dic.keys():
                print('----')
                try:
                    print(self.__default_property.is_leaf)
                    print(child)
                    print(default_dic[child])
                    print(config_dic[child])
                except:
                    pass

                if config_dic != None:
                    print('!-!-!')
                    print(config_dic)
                    print('####')
                    pass

                    try:
                        config_dic[child]
                        #self.children[child] = PropertyBase(default_dic[child], config_dic[child], child)
                        tmp[child] = PropertyBase(default_dic[child], config_dic[child], child)
                    except (AttributeError, TypeError, KeyError):
                        #self.children[child] = PropertyBase(default_dic[child], config_dic, child)
                        tmp[child] = PropertyBase(default_dic[child], config_dic, child)

            self.children[default_section] = tmp
"""

"""
        def PF(top_v):
            tmp = {}
            default_property = DefaultParser(top_v)
            print(top_v)
            if not default_property.is_leaf:
                tmp['branch_properties'] = default_property
                print(top_v)
                for k, v in top_v.items():
                    print("key is %s" % str(k))
                    tmp[k] = PF(v)
                return tmp
            else:
                return default_property
"""

"""
    def get_config(self):
        return self.conf_v

    def get_default(self):




    def __finditem(self, obj, key):
        if key in obj: return obj[key]
        for k, v in obj.items():
            if isinstance(v,dict):
                item = _finditem(v, key)
                if item is not None:
                    return item
"""
