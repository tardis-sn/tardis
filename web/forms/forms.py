#TODO: quantity_range_sampled
#TODO: Add support for container-property

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
    'string': fields.TextField,
    'list': fields.TextField,
    'bool': fields.BooleanField,
}

#------------------------------------------------------------------------------

def createField(schema,parent=None):
    description = {}
    validation = []
    description['help_text'] = schema['help']
    if 'default' in schema and schema['default']!='None':
        description['default'] = schema['default']
    else:
        description['default'] = ''
    if parent:
        description['parent']=parent
    if 'mandatory' in schema and schema['mandatory']:
        validation = [validators.DataRequired()]
    if 'file' in schema and schema['file']:
        form_field = fields.FileField
    elif 'allowed_value' in schema:
        form_field = fields.SelectField
        choices = [(i,i) for i in schema['allowed_value'].split()]
        return form_field(choices=choices,validators=validation,description=description)
    else:
        form_field = field_type[schema['property_type']]
    return form_field(validators=validation,description=description)

def autofill_fields(yml_field,single_item):
    def run(cls):
        print yml_field
        print
        if single_item:
            setattr(cls, single_item, createField(yml_field))
        else:
            for i in yml_field:
                if 'property_type' in yml_field[i]:
                    setattr(cls, i, createField(yml_field[i]))
                else:
                    setattr(cls, i, fields.BooleanField())
                    for j in yml_field[i]:
                        setattr(cls, j, createField(yml_field[i][j],parent=i))

        return cls
    return run

#------------------------------------------------------------------------------

@autofill_fields(config_data['supernova'],single_item=False)
class SupernoveForm(Form):
    pass

@autofill_fields(config_data['atom_data'],single_item='atom_data')
class AtomForm(Form):
    pass

@autofill_fields(config_data['plasma'],single_item=False)
class PlasmaForm(Form):
    pass

@autofill_fields(config_data['spectrum'],single_item='spectrum')
class SpectrumForm(Form):
    pass

class AbundanceForm(Form):
    pass

class StructureForm(Form):
    pass

class MonteCarloForm(Form):
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
