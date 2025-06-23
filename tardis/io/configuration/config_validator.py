"""
Configuration validation module for TARDIS.

This module provides functionality to validate TARDIS configuration files
against JSON schemas with support for cross-schema references and custom
type checking for Astropy Quantity objects.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

import yaml
from astropy.units.quantity import Quantity
from jsonschema import validators
from jsonschema.validators import Draft7Validator
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT7

from tardis.io.util import YAMLLoader

if TYPE_CHECKING:
    from typing import Any

CONFIGURATION_DIR = Path(__file__).resolve().parent
SCHEMA_DIR = CONFIGURATION_DIR / "schemas"
CONFIG_SCHEMA_FNAME = SCHEMA_DIR / "base.yml"


def extend_with_default(
    validator_class: type[Draft7Validator],
) -> type[Draft7Validator]:
    """
    Extend a jsonschema validator to automatically set default values on properties.

    By default, jsonschema ignores default values defined in schemas. This function
    creates an extended validator that will set default values during validation.

    Parameters
    ----------
    validator_class : type[Draft7Validator]
        The jsonschema validator class to extend (e.g., Draft7Validator)

    Returns
    -------
    type[Draft7Validator]
        Extended validator class that sets default values during validation

    Notes
    -----
    The extended validator also properly handles default values in schemas that
    use the oneOf keyword and validates that default values are of correct types.
    """
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(
        validator: Any,
        properties: dict[str, Any],
        instance: dict[str, Any],
        schema: dict[str, Any],
    ) -> Any:
        """Set default values for properties that don't exist in the instance."""
        # This validator also checks if default values are of the correct type
        # and properly sets default values on schemas that use the oneOf keyword
        if not list(validate_properties(validator, properties, instance, schema)):
            for property_name, subschema in properties.items():
                if "default" in subschema:
                    instance.setdefault(property_name, subschema["default"])

        yield from validate_properties(validator, properties, instance, schema)

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


DefaultDraft7Validator = extend_with_default(Draft7Validator)


def _create_schema_registry() -> Registry:
    """
    Create a registry containing all schema files for reference resolution.

    Loads all YAML schema files from the schemas directory and creates a
    referencing Registry that can be used to resolve $ref references between
    schema files during validation.

    Returns
    -------
    Registry
        A referencing Registry containing all schema files with their
        filenames as URIs for cross-reference resolution

    Notes
    -----
    Schema files are loaded using the YAMLLoader and registered with
    Draft 7 JSON Schema specification for proper validation.
    """
    registry = Registry()

    # Load all schema files in the schemas directory
    for schema_file in SCHEMA_DIR.glob("*.yml"):
        with open(schema_file) as f:
            schema_content = yaml.load(f, Loader=YAMLLoader)

        # Create a resource with the filename as the URI, explicitly specifying Draft 7
        resource = Resource.from_contents(schema_content, default_specification=DRAFT7)
        registry = registry.with_resource(uri=schema_file.name, resource=resource)

    return registry


def is_quantity(checker: Any, instance: Any) -> bool:
    """
    Check if the provided instance is of type astropy.units.quantity.Quantity.

    This function is used as a custom type checker for jsonschema validation
    to properly handle Astropy Quantity objects in configuration validation.

    Parameters
    ----------
    checker : Any
        Object of `TypeChecker`. Passed by jsonschema internally.
    instance : Any
        The instance to be checked for Quantity type.

    Returns
    -------
    bool
        True if the instance is of type astropy.units.quantity.Quantity,
        False otherwise.
    """
    return isinstance(instance, Quantity)


def validate_dict(
    config_dict: dict[str, Any],
    schemapath: Path = CONFIG_SCHEMA_FNAME,
    validator: type[Draft7Validator] = DefaultDraft7Validator,
) -> dict[str, Any]:
    """
    Validate a configuration dictionary against a JSON schema.

    This function validates a configuration dictionary using a JSON schema,
    with support for cross-file schema references and custom type checking
    for Astropy Quantity objects. Default values are automatically set
    according to the schema definition.

    Parameters
    ----------
    config_dict : dict[str, Any]
        The configuration dictionary to validate.
    schemapath : Path, optional
        Path to the main schema file. Defaults to CONFIG_SCHEMA_FNAME.
    validator : type[Draft7Validator], optional
        The validator class to use. Defaults to DefaultDraft7Validator.

    Returns
    -------
    dict[str, Any]
        A validated and potentially modified copy of the input dictionary
        with default values set according to the schema.

    Raises
    ------
    jsonschema.ValidationError
        If the configuration dictionary does not conform to the schema.
    jsonschema.SchemaError
        If the schema itself is invalid.

    Notes
    -----
    This function creates a deep copy of the input dictionary before
    validation to avoid modifying the original data.
    """
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
    configpath: Path,
    schemapath: Path = CONFIG_SCHEMA_FNAME,
    validator: type[Draft7Validator] = DefaultDraft7Validator,
) -> dict[str, Any]:
    """
    Validate a YAML configuration file against a JSON schema.

    This function loads a YAML configuration file and validates it using
    the validate_dict function. It provides a convenient interface for
    validating configuration files directly from disk.

    Parameters
    ----------
    configpath : Path
        Path to the YAML configuration file to validate.
    schemapath : Path, optional
        Path to the main schema file. Defaults to CONFIG_SCHEMA_FNAME.
    validator : type[Draft7Validator], optional
        The validator class to use. Defaults to DefaultDraft7Validator.

    Returns
    -------
    dict[str, Any]
        A validated configuration dictionary with default values set
        according to the schema.

    Raises
    ------
    FileNotFoundError
        If the configuration file cannot be found.
    yaml.YAMLError
        If the YAML file cannot be parsed.
    jsonschema.ValidationError
        If the configuration does not conform to the schema.
    jsonschema.SchemaError
        If the schema itself is invalid.

    Notes
    -----
    The YAML file is loaded using the custom YAMLLoader which may have
    specific behavior for handling certain YAML constructs.
    """
    with open(configpath) as f:
        config = yaml.load(f, Loader=YAMLLoader)
    return validate_dict(config, schemapath, validator)
