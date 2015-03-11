#TODO: Added a property for label names in the YAML File, to make the web-app more clear.

from wtforms import validators, fields
from wtforms_tornado import Form
import yaml
import json
from collections import OrderedDict

mandatory = [validators.DataRequired()]

#------------------------------------------------------------------------------

def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)

config_data = ordered_load(file('config.yml', 'r'))

field_type = {
    'quantity': fields.TextField,
    'quantity_range_sampled': fields.TextField,
    'float': fields.TextField,
    'int': fields.TextField,
    'string': fields.TextField,
    'list': fields.TextField,
    'bool': fields.BooleanField,
    'container-declaration': fields.SelectField
}

association_dict = {}

#------------------------------------------------------------------------------

def createField(schema,parent='',name=None):
    description = {}
    validation = []
    if 'help' in schema:
        description['help_text'] = schema['help']
    if 'default' in schema and schema['default']!='None':
        description['default'] = schema['default']
    else:
        description['default'] = ''
    if name in association_dict and parent==association_dict[name]['parent']:
        try:
            description['association']=association_dict[name]['association']
        except:
            pass

    if parent:
        description['parent']=parent
    if 'mandatory' in schema and schema['mandatory']:
        validation = [validators.DataRequired()]
    # if 'file' in schema and schema['file']:
        # form_field = fields.FileField
    if 'allowed_value' in schema:
        form_field = fields.SelectField
        choices = [(i,i) for i in schema['allowed_value']]
        return form_field(choices=choices,validators=validation,description=description)
    else:
        form_field = field_type[schema['property_type']]
        if schema['property_type']=='container-declaration':
            choices=[(i,i) for i in schema['containers']]
            children = [i  for i in schema if (i[0]=='_' or i[0]=='+')]
            for i in children:
                for j in schema[i]:
                    if j in association_dict:
                        association_dict[j]['association'] = association_dict[j]['association']+" "+i[1:]
                    else:
                        association_dict[j] = {}
                        if parent:
                            association_dict[j]['parent'] = parent
                            association_dict[j]['association'] = i[1:] + " " +parent+'|'+ name
                        else:
                            association_dict[j]['parent'] = parent
                            association_dict[j]['association'] = i[1:] + " " +name

            return form_field(choices=choices,validators=validation,description=description)
    return form_field(validators=validation,description=description)

def populate_fields(yml_field,single_item,parent=''):
    def run(cls):
        global association_dict
        if single_item:
            setattr(cls, single_item, createField(yml_field))
        else:
            for i in yml_field:
                try:
                    association = association_dict[i]['association']
                except:
                    association = ''
                if i=='property_type' and yml_field[i]=='container-property':
                    pass
                elif 'property_type' in yml_field[i]:

                    if 'property_type' in yml_field[i] and yml_field[i]['property_type']=='container-property':
                        if association:
                            setattr(cls, i, fields.BooleanField(description={'association':association}))
                        else:
                            setattr(cls, i, fields.BooleanField())
                        populate_fields(yml_field[i],single_item,parent=i)(cls)
                        association_dict = {}
                    else:
                        if parent:
                            setattr(cls, parent+'|'+i, createField(yml_field[i],parent=parent,name=i))
                        else:
                            setattr(cls, i, createField(yml_field[i],parent=parent,name=i))
                else:
                    if parent:
                        setattr(cls, i, fields.BooleanField(description={'association':association,'parent':parent}))
                    else:
                        setattr(cls, i, fields.BooleanField(description={'association':association}))
                    populate_fields(yml_field[i],single_item,parent=i)(cls)
        return cls
    return run

#------------------------------------------------------------------------------

@populate_fields(config_data['supernova'],single_item=False)
class SupernoveForm(Form):
    pass

@populate_fields(config_data['atom_data'],single_item='atom_data')
class AtomForm(Form):
    pass

@populate_fields(config_data['plasma'],single_item=False)
class PlasmaForm(Form):
    pass

@populate_fields(config_data['model']['abundances'],single_item=False)
class AbundanceForm(Form):
    uniform_abundances = fields.TextAreaField(description={'help_text':'Insert Uniform abundances of all the shells, in the format: C: 0.01 O: 0.01 etc...'})

@populate_fields(config_data['model']['structure'],single_item=False)
class StructureForm(Form):
    pass

@populate_fields(config_data['montecarlo'],single_item=False)
class MonteCarloForm(Form):
    pass

@populate_fields(config_data['spectrum'],single_item='spectrum')
class SpectrumForm(Form):
    pass
#------------------------------------------------------------------------------

class TardisForm(Form):
    supernova = fields.FormField(SupernoveForm)
    atom_data = fields.FormField(AtomForm)
    plasma = fields.FormField(PlasmaForm)
    structure_model = fields.FormField(StructureForm)
    abundances_model = fields.FormField(AbundanceForm)
    montecarlo = fields.FormField(MonteCarloForm)
    spectrum = fields.FormField(SpectrumForm)
