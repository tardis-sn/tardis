from copy import deepcopy
from pathlib import Path

import yaml
from astropy.units.quantity import Quantity
from jsonschema import Draft4Validator, RefResolver, validators

from tardis.io.util import YAMLLoader

CONFIGURATION_DIR = Path(__file__).resolve().parent
SCHEMA_DIR = CONFIGURATION_DIR / "schemas"
CONFIG_SCHEMA_FNAME = SCHEMA_DIR / "base.yml"


def extend_with_default(validator_class):
    """
    Extend a `jsonschema.IValidator` to also set default
    values on properties. By default jsonschema ignores
    default values.

    Parameters
    ----------
    validator_class :
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
        raise Exception(f"Not a file URL: {path}")
    with open(path[len("file://") :]) as f:
        return yaml.load(f, Loader=YAMLLoader)


def is_quantity(checker, instance):
    """
    Check if the provided instance is of type astropy.units.quantity.Quantity

    Parameters
    ----------
    checker:
        Object of `TypeChecker`. Passed by jsonschema internally.

    instance:
        The instance to be checked.

    Returns
    -------
    bool: True if the instance is of type astropy.units.quantity.Quantity else False
    """
    return isinstance(instance, Quantity)


def validate_dict(
    config_dict,
    schemapath=CONFIG_SCHEMA_FNAME,
    validator=DefaultDraft4Validator,
):
    with open(schemapath) as f:
        schema = yaml.load(f, Loader=YAMLLoader)
    schemaurl = f"file://{schemapath}"
    handlers = {"file": _yaml_handler}
    resolver = RefResolver(schemaurl, schema, handlers=handlers)
    validated_dict = deepcopy(config_dict)
    custom_type_checker = validator.TYPE_CHECKER.redefine(
        "quantity", is_quantity
    )
    custom_validator = validators.extend(
        validator, type_checker=custom_type_checker
    )
    custom_validator(
        schema=schema,
        resolver=resolver,
    ).validate(validated_dict)
    return validated_dict


def validate_yaml(
    configpath, schemapath=CONFIG_SCHEMA_FNAME, validator=DefaultDraft4Validator
):
    with open(configpath) as f:
        config = yaml.load(f, Loader=YAMLLoader)
    return validate_dict(config, schemapath, validator)
