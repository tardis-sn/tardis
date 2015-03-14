import urllib2
import yaml

from django import forms
from django.utils.encoding import smart_text

_URL = 'https://raw.githubusercontent.com/tardis-sn/tardis/master/tardis/data/tardis_config_definition.yml'

_FIELD_MAP = {
    'string': forms.CharField,
    'list': forms.CharField,
    'quantity_range_sampled': forms.CharField,
    'container-property': forms.CharField,
    'quantity': forms.FloatField,
    'bool': forms.BooleanField,
    'float': forms.FloatField,
    'int': forms.IntegerField
}


def generate_schema():
    """
    ```generate_schema``` generates a dictionary
    from the the configuration scheme.
    The configuration schema is retrived from _URL.
    """
    return substitute_fields(yaml.load(open('tardis/tardis_config_definition.yml')))


def substitute_fields(field_definition):
    """
    """
    for (field_name, field) in field_definition.iteritems():
        dtype = field.get('property_type', None)
        if dtype is None or dtype == 'container-property':
            substitute_fields(field)
        else:
            if 'default' in field:
                if field['default'] == 'None':
                    del field['default']
            if dtype in _FIELD_MAP:
                field_definition[field_name] = _FIELD_MAP[dtype](
                    # label=smart_text(field_name),
                    initial=field.get('default', None),
                    required=field.get('mandatory', None),
                    help_text=field.get('help', None)
                )
    return field_definition
