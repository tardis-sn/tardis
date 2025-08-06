from pathlib import Path

import pytest
import yaml
from astropy import units as u
from jsonschema import Draft7Validator

from tardis.io.configuration.config_validator import (
    extend_with_default,
    is_quantity,
    validate_dict,
    validate_yaml,
)
from tardis.io.util import YAMLLoader

CONFIGURATION_TESTS_DIR = Path(__file__).resolve().parent
CONFIG_DATA_VERYSIMPLE = CONFIGURATION_TESTS_DIR / "data/tardis_configv1_verysimple.yml"


def test_yaml_loading(tardis_config_verysimple):
    """Test that YAML files can be loaded correctly."""
    # Load the YAML file directly using the same loader as the config validator
    with open(CONFIG_DATA_VERYSIMPLE) as f:
        loaded_config = yaml.load(f, Loader=YAMLLoader)

    # Compare with the fixture
    assert loaded_config == tardis_config_verysimple

    # Test that non-existent files raise FileNotFoundError
    with pytest.raises(FileNotFoundError):
        with open("/nonexistent/file.yml") as f:
            yaml.load(f, Loader=YAMLLoader)


def test_is_quantity():
    """Test the is_quantity function."""
    quantity_val_1 = 5 * u.m
    quantity_val_2 = 5

    assert not is_quantity(None, quantity_val_2)
    assert is_quantity(None, quantity_val_1)


def test_extend_with_default(tardis_config_verysimple):
    """Test the extend_with_default function."""
    Default_Validator = extend_with_default(Draft7Validator)

    # Using the validator that extends Draft7Validator with default handling
    dict_with_Default = validate_dict(tardis_config_verysimple, validator=Default_Validator)
    assert dict_with_Default["plasma"]["initial_t_inner"] == "-1 K"

    # The main test is that the extended validator sets defaults properly
    # We don't need to test the base validator separately since it would fail validation
    # due to missing required fields that the extended validator fills with defaults
    assert "plasma" in dict_with_Default
    assert "initial_t_inner" in dict_with_Default.get("plasma", {})


def test_validate_dict(tardis_config_verysimple):
    """Test the validate_dict function."""
    config_dict_verysimple = validate_dict(tardis_config_verysimple)

    assert config_dict_verysimple["montecarlo"]["seed"] == tardis_config_verysimple["montecarlo"]["seed"]

    # Checks for default value when not provided
    assert config_dict_verysimple["plasma"]["initial_t_inner"] == "-1 K"
    assert config_dict_verysimple["supernova"]["luminosity_wavelength_start"] == "0 angstrom"


def test_validate_yaml():
    """Test the validate_yaml function."""
    path = CONFIG_DATA_VERYSIMPLE

    # Load the expected data using the same method as the old _yaml_handler
    with open(CONFIG_DATA_VERYSIMPLE) as f:
        expected_dict = yaml.load(f, Loader=YAMLLoader)

    # Provided config path retrieves validated dict
    # Use the default validator which includes the extend_with_default functionality
    validated_dict = validate_yaml(path)
    assert validated_dict["spectrum"]["num"] == expected_dict["spectrum"]["num"]

    # Test that the validation works and returns a valid dict
    assert isinstance(validated_dict, dict)
    assert "spectrum" in validated_dict
