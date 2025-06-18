from copy import deepcopy
from pathlib import Path

import yaml
from astropy.units.quantity import Quantity
from jsonschema import Draft7Validator, validators
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT7

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

        yield from validate_properties(
            validator,
            properties,
            instance,
            schema,
        )

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


DefaultDraft7Validator = extend_with_default(Draft7Validator)


def _create_schema_registry():
    """Create a registry containing all schema files for reference resolution."""
    registry = Registry()

    # Load all schema files in the schemas directory
    for schema_file in SCHEMA_DIR.glob("*.yml"):
        with open(schema_file) as f:
            schema_content = yaml.load(f, Loader=YAMLLoader)

        # Create a resource with the filename as the URI, explicitly specifying Draft 7
        resource = Resource.from_contents(schema_content, default_specification=DRAFT7)
        registry = registry.with_resource(uri=schema_file.name, resource=resource)

    return registry


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
    validator=DefaultDraft7Validator,
):
    with open(schemapath) as f:
        schema = yaml.load(f, Loader=YAMLLoader)

    # Create registry for schema references
    registry = _create_schema_registry()

    validated_dict = deepcopy(config_dict)
    custom_type_checker = validator.TYPE_CHECKER.redefine(
        "quantity", is_quantity
    )
    custom_validator = validators.extend(
        validator, type_checker=custom_type_checker
    )

    # Create validator with registry for reference resolution
    validator_instance = custom_validator(schema=schema, registry=registry)
    validator_instance.validate(validated_dict)
    return validated_dict


def validate_yaml(
    configpath, schemapath=CONFIG_SCHEMA_FNAME, validator=DefaultDraft7Validator
):
    with open(configpath) as f:
        config = yaml.load(f, Loader=YAMLLoader)
    return validate_dict(config, schemapath, validator)
