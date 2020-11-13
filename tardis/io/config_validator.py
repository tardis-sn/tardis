import os
import yaml
from copy import deepcopy
from jsonschema import Draft4Validator, RefResolver, validators
from astropy.units.quantity import Quantity
from tardis.io.util import YAMLLoader

base_dir = os.path.abspath(os.path.dirname(__file__))
schema_dir = os.path.join(base_dir, "schemas")
config_schema_file = os.path.join(schema_dir, "base.yml")


def extend_with_default(validator_class):
    """
    Extend a `jsonschema.IValidator` to also set default
    values on properties. By default jsonschema ignores
    default values.

    Parameters
    ----------
    validator_class:
        The `jsonschema.IValidator` class to extend

    Returns
    -------
    The extended `jsonschema.IValidator`

    """
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        # This validator also checks if default values
        # are of the correct type and properly sets default
        # values on schemas that use the oneOf keyword
        if not list(
            validate_properties(validator, properties, instance, schema)
        ):
            for property, subschema in properties.items():
                if "default" in subschema:
                    instance.setdefault(property, subschema["default"])

        for error in validate_properties(
            validator,
            properties,
            instance,
            schema,
        ):
            yield error

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


DefaultDraft4Validator = extend_with_default(Draft4Validator)


def _yaml_handler(path):
    if not path.startswith("file://"):
        raise Exception("Not a file URL: {}".format(path))
    with open(path[len("file://") :]) as f:
        return yaml.load(f, Loader=YAMLLoader)


def validate_dict(
    config_dict, schemapath=config_schema_file, validator=DefaultDraft4Validator
):
    with open(schemapath) as f:
        schema = yaml.load(f, Loader=YAMLLoader)
    schemaurl = "file://" + schemapath
    handlers = {"file": _yaml_handler}
    resolver = RefResolver(schemaurl, schema, handlers=handlers)
    validated_dict = deepcopy(config_dict)
    validator(
        schema=schema, types={"quantity": (Quantity,)}, resolver=resolver
    ).validate(validated_dict)
    return validated_dict


def validate_yaml(
    configpath, schemapath=config_schema_file, validator=DefaultDraft4Validator
):
    with open(configpath) as f:
        config = yaml.load(f, Loader=YAMLLoader)
    return validate_dict(config, schemapath, validator)
