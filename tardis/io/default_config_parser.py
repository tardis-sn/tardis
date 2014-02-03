import yaml as yaml
from astropy import units
import re

class DefaultParser:
    """Not invented here syndrome"""

    __check = {}
    __convert = {}
    __list_of_leaf_types = []


    def __init__(self, default_dict):
        
        self.__register_leaf('list')
        self.__register_leaf('int')
        self.__register_leaf('float')
        self.__register_leaf('quantity')
        self.__register_leaf('string')
        self.__register_leaf('container-declaration')
        self.__mandatory = False


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
            self.__allowed_value = self.__convert_av_in_pt(default_dict['allowed_value'], self.__property_type)

        if 'allowed_type' in default_dict:
            self.__allowed_type = default_dict['allowed_type']
            self.__lower, self.__upper = self.__parse_allowed_type(self.__allowed_type)

        if 'default' in default_dict:
            self.set_default(default_dict['default'])

        if 'mandatory' in default_dict:
            self.__mandatory =  default_dict['mandatory']

    

        self.is_leaf = self.__is_leaf(self.__property_type)

    def get_default(self):
        return self.__default_value

    def set_default(self, value):
        if value != None:
            if self.is_valid(value):
                self.__default_value = value
            else:
                raise ValueError('Default value violates property constraint.')
        else:
            self.__default_value = None

    def is_mandatory(self):
        return self.__mandatory

    def has_default(self):
        try:
            if self.__default_value:
                return True
            else:
                return False
        except NameError:
            pass

    def set_path_in_dic(self, path):
        self.__path = path
        
    def get_path_in_dict(self):
        return self.__path

    def set_config_value(self, value):
        self.__config_value = value

    def get_value(self):
        
        if self.__config_value is not None and  self.is_valid(self.__config_value):
            return self.__config_value
        else:
            if self.has_default():
               return self.__default_value
            else:
               raise ValueError('No default value given.')

    def is_container(self):
        return self.__is_container(None)

    def get_container_dic(self):
        if self.__is_container(None):
            return self.__container_dic

    def update_container_dic(self, container_dic, current_entry_name):
        if reduce(lambda a,b: a or b, [container_dic.has_key(i) for i in ['and','or']], True):
            if 'or' in container_dic:
                if current_entry_name in container_dic['or']:
                    container_dic['or'] = []
                    return container_dic
            if 'and' in container_dic:
                if current_entry_name in container_dic['and']:
                    current_entry_name['and'].remove(current_entry_name)
                    return container_dic






    def is_valid(self, value):
        if not self.__check[self.__property_type](self,value):
            return False
        if self.__allowed_value:
            if not self.__is_allowed_value(self,value, self.__allowed_value):
                return False
        if self.__allowed_type:
            if not self.__check_value(self,value, self.__lower, self.__upper):
                return False
        return True

    def __register_leaf(self,type_name):
        print(type_name)
        if not type_name in self.__list_of_leaf_types:
            self.__list_of_leaf_types.append(type_name)

    def __is_leaf(self, type_name):
        return type_name in self.__list_of_leaf_types

    def __is_container(self, value):
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
            

    def __convert_av_in_pt(self,allowed_value, property_type):
        """
        Converts the allowed values to the property type.
        """
        if not len([]) == 0:
            return [self.__convert[property_type](a) for a in property_value]
        else:
            return []


    def __is_allowed_value(self, value,  allowed_value):
        if value in _allowed_value:
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
                        tmp[child]= PropertyBase(default_dic[child], config_dic[child], child)
                    except (AttributeError,TypeError, KeyError):
                        #self.children[child] = PropertyBase(default_dic[child], config_dic, child)
                        tmp[child]= PropertyBase(default_dic[child], config_dic, child)
                    
            self.children[default_section] = tmp

class Container(DefaultParser):

    def __init__(self, container_default_dict, container_dict):
        #self.__register_leaf('list')
        #self.__register_leaf('int')
        #self.__register_leaf('float')
        #self.__register_leaf('quantity')
        #self.__register_leaf('string')

        self.__allowed_value = None
        self.__allowed_type = None
        self.__config_value = None
        self.__path = None

        self.__property_type =  'container-property'

        self.__default_container = {}
        self.__config_container = {} 

        #check if it is a valid default container
        if not 'type' in container_default_dict:
            raise ValueError('The given default contaienr is no valid')

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
            self.__selected_container =  container_dict['type']
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
            print('START NEW PARSE')
            tmp_conf_ob = {}
            path_in_dic = []
            tmp_conf_val = {}
            if isinstance(top_default,dict):
                print(top_default)
                default_property = DefaultParser(top_default)
                print(default_property.is_leaf)
                if not default_property.is_leaf:
                    print(top_default.items())
                    for k,v in top_default.items():
                        print('>>--<<')
                        print(top_config)
                        print(k)
                        print(v)
                        tmp_conf_ob[k], tmp_conf_val[k]  = parse_container_items(v, top_config, k, path + [k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
                        print('>>>>>>>>>>>')
                        print(path)
                        print(top_config)
                        conf_value = get_property_by_path(top_config, path)
                        print('conf_value: %s'%conf_value)
                    except:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)

                    return default_property, default_property.get_value()

        def get_property_by_path(conf_dict, path):
            for key in path:
                conf_dict = conf_dict[key]
            return conf_dict
        
        for item in necessary_items:
            if not item in container_dict.keys():
                raise ValueError('Entry %s is missing in container'%str(item))
            self.__default_container, self.__config_container = parse_container_items(container_default_dict[item], container_dict[item], item, [])
            pass
                #go through  all items and create an conf object thereby check the conf

        self.__container_ob = self.__default_container
        self.__conf = self.__config_container

    def get_container_ob(self):
        return self.__container_ob

    def get_container_conf(self):
        return self.__conf




class Config:
    
    def __init__(self, default_conf, conf):
        self.mandatories = {}
        self.fulfilled = {}
        self.__create_default_conf(default_conf)
        self.__parse_config(default_conf, conf)

    def __mandatory_key(self, path):
        return ':'.join(path)

    def register_mandatory(self, name, path):
        self.mandatories[self.__mandatory_key(path)] = name

    def deregister_mandatory(self, name, path):
        self.fulfilled[self.__mandatory_key(path)] = name

    def is_mandatory_fulfilled(self):
        if len(set(self.mandatories.keys()) - set(self.fulfilled.keys()))<=0:
            return True
        else:
            return False





    def __parse_config(self, default_conf, conf,):

        def PF( top_v):
            tmp = {}
            default_property = DefaultParser(top_v)
            print(top_v)
            if not default_property.is_leaf:
                tmp['branch_properties'] = default_property
                print(top_v)
                for k,v in top_v.items():
                    print("key is %s"%str(k))
                    tmp[k] = PF(v)
                return tmp
            else:
                return default_property

        def finditem( obj, key):
            if key in obj: return obj[key]
            for k, v in obj.items():
                if isinstance(v,dict):
                    item = finditem(v, key)
                    if item is not None:
                        return item

        def is_path_valid(conf_dict, path):
            try:
                for key in path:
                    conf_dict = conf_dict[key]
                return True
            except KeyError:
                return False


        def get_property_by_path(conf_dict, path):
            for key in path:
                conf_dict = conf_dict[key]
            return conf_dict
        def recursive_parser(top_v, conf, level_name, path):
            tmp_conf_ob = {}
            path_in_dic = []
            tmp_conf_val = {}
            if isinstance(top_v,dict):
                default_property = DefaultParser(top_v)
                if default_property.is_mandatory():
                    self.register_mandatory(self, path)
                self.deregister_mandatory(self, path)

                if default_property.is_container():
                    container_conf = get_property_by_path(conf, path)
                    ccontainer = Container(top_v, container_conf)
                    return ccontainer.get_container_ob(), ccontainer.get_container_conf()
                elif not default_property.is_leaf:
                    for k,v in top_v.items():
                        tmp_conf_ob[k], tmp_conf_val[k]  = recursive_parser(v, conf, k, path + [k])
                    return tmp_conf_ob, tmp_conf_val
                else:
                    default_property.set_path_in_dic(path)
                    try:
                        conf_value = get_property_by_path(conf, path)
                    except:
                        conf_value = None

                    if conf_value is not None:
                        default_property.set_config_value(conf_value)
                    return default_property, default_property.get_value()

        
        self.__conf_o, self.__conf_v = recursive_parser(default_conf, conf, 'main', [])

    def __create_default_conf(self, default_conf):
        def recursive_default_parser(top_v,level_name, path):
            tmp_default = {}
            path_in_dic = []
            if isinstance(top_v,dict):
                default_property = DefaultParser(top_v)
                if not default_property.is_container():
                    if not default_property.is_leaf:
                        for k,v in top_v.items():
                            tmp_default[k] = recursive_default_parser(v, k, path + [k])
                        return tmp_default
                    else:
                        default_property.set_path_in_dic(path)
                        if default_property.has_default():
                            return default_property.get_default()
                        else:
                            return None
        self.__default_config = recursive_default_parser(default_conf, 'main', [])
        
    def get_config(self):
        return self.__conf_v

    def get_default_config(self):
        return self.__default_config

    def get_config_object(self):
        return self.__conf_o

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
