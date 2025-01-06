import pytest 
from jsonschema import Draft4Validator
from pathlib import Path
from astropy import units as u
from jsonschema.exceptions import ValidationError

from tardis.io.configuration.config_validator import (
    _yaml_handler, 
    validate_dict, 
    validate_yaml, 
    extend_with_default,
    is_quantity
)

CONFIGURATION_TESTS_DIR = Path(__file__).resolve().parent
CONFIG_DATA_VERYSIMPLE = CONFIGURATION_TESTS_DIR / "data/tardis_configv1_verysimple.yml"


def test_yaml_handler(tardis_config_verysimple):
    res = _yaml_handler(f"file://{CONFIG_DATA_VERYSIMPLE}")
    assert res == tardis_config_verysimple

    with pytest.raises(Exception):
        _yaml_handler("/mock/example.yml")

    with pytest.raises(FileNotFoundError):
        _yaml_handler("file:///example.yml")


def test_is_quantity():
    quantity_val_1 = 5 * u.m
    quantity_val_2 = 5

    assert is_quantity(None, quantity_val_2) == False
    assert is_quantity(None, quantity_val_1) == True


def test_extend_with_default(tardis_config_verysimple):
    Default_Validator = extend_with_default(Draft4Validator)

    # Using the validator that extended from Draft4
    dict_with_Default = validate_dict(tardis_config_verysimple, validator=Default_Validator)
    assert dict_with_Default["plasma"]["initial_t_inner"] == "-1 K"

    with pytest.raises(ValidationError):
        dict_with_Draft4 = validate_dict(tardis_config_verysimple, validator=Draft4Validator)


def test_validate_dict(tardis_config_verysimple):
    config_dict_verysimple = validate_dict(tardis_config_verysimple)

    assert config_dict_verysimple["montecarlo"]["seed"] == tardis_config_verysimple["montecarlo"]["seed"]

    # Checks for default value when not provided
    assert config_dict_verysimple["plasma"]["initial_t_inner"] == "-1 K"
    assert config_dict_verysimple["supernova"]["luminosity_wavelength_start"] == "0 angstrom"

    # ValidationError because Draft4Validator cannot process default values
    # Context: cannot assign default value for model_isotope_time_0 in abundances which is a required field
    with pytest.raises(ValidationError):
        validate_dict(tardis_config_verysimple, validator=Draft4Validator)


def test_validate_yaml():
    path = CONFIG_DATA_VERYSIMPLE
    dict = _yaml_handler(f"file://{CONFIG_DATA_VERYSIMPLE}")

    # Provided config path retrieves validated dict
    assert validate_yaml(path)["spectrum"]["num"] == dict["spectrum"]["num"] 

    with pytest.raises(ValidationError):
        validate_yaml(path, validator=Draft4Validator)
