# tests for the config reader module
import os
from attr import validate
import pytest
import pandas as pd
from numpy.testing import assert_almost_equal
from jsonschema.exceptions import ValidationError

from tardis.io import config_reader
from tardis.io.config_reader import Configuration


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.abspath(os.path.join(data_dir, "data", filename))


def test_convergence_section_parser():
    test_convergence_section = {
        "type": "damped",
        "lock_t_inner_cyles": 1,
        "t_inner_update_exponent": -0.5,
        "damping_constant": 0.5,
        "threshold": 0.05,
        "fraction": 0.8,
        "hold_iterations": 3,
        "t_rad": {"damping_constant": 1.0},
    }

    parsed_convergence_section = config_reader.parse_convergence_section(
        test_convergence_section
    )

    assert_almost_equal(
        parsed_convergence_section["t_rad"]["damping_constant"], 1.0
    )

    assert_almost_equal(
        parsed_convergence_section["w"]["damping_constant"], 0.5
    )


def test_from_config_dict(tardis_config_verysimple):
    conf = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )
    assert conf.config_dirname == "test"

    assert_almost_equal(
        conf.spectrum.start.value,
        tardis_config_verysimple["spectrum"]["start"].value,
    )

    assert_almost_equal(
        conf.spectrum.stop.value,
        tardis_config_verysimple["spectrum"]["stop"].value,
    )

    tardis_config_verysimple["spectrum"]["start"] = "Invalid"
    with pytest.raises(ValidationError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )


def test_config_hdf(hdf_file_path, tardis_config_verysimple):
    expected = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )
    expected.to_hdf(hdf_file_path, overwrite=True)
    actual = pd.read_hdf(hdf_file_path, key="/simulation/config")
    expected = expected.get_properties()["config"]
    assert actual[0] == expected[0]


def test_model_section_config(tardis_config_verysimple):
    """
    Configuration Validation Test for Model Section of the Tardis Config YAML File

    Validates:
        Density: branch85_w7
        Velocity (Start < End)

    Parameter
    ---------
        `tardis_config_verysimple` : YAML File

    Result
    ------
        Assertion based on validation for specified values
    """
    conf = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )

    assert conf.model.structure.density.type == "branch85_w7"

    tardis_config_verysimple["model"]["structure"]["velocity"][
        "start"
    ] = "2.0e4 km/s"
    tardis_config_verysimple["model"]["structure"]["velocity"][
        "stop"
    ] = "1.1e4 km/s"
    with pytest.raises(ValueError) as ve:
        if (
            conf.model.structure.velocity.start
            < conf.model.structure.velocity.stop
        ):
            raise ValueError("Stop Value must be greater than Start Value")
    assert ve.type is ValueError


def test_supernova_section_config(tardis_config_verysimple):
    """
    Configuration Validation Test for Supernova Section of the Tardis Config YAML File

    Validates:
        Time of Explosion (Must always be positive)
        Luminosity Wavelength Limits (Start < End)

    Parameter
    ---------
        `tardis_config_verysimple` : YAML File

    Result
    ------
        Assertion based on validation for specified values
    """
    conf = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )
    tardis_config_verysimple["supernova"]["time_explosion"] = "-10 day"
    tardis_config_verysimple["supernova"][
        "luminosity_wavelength_start"
    ] = "15 angstrom"
    tardis_config_verysimple["supernova"][
        "luminosity_wavelength_end"
    ] = "0 angstrom"
    with pytest.raises(ValueError) as ve:
        if conf.supernova.time_explosion.value > 0:
            raise ValueError("Time of Explosion cannot be negative")
    assert ve.type is ValueError

    with pytest.raises(ValueError) as ve:
        if (
            conf.supernova.luminosity_wavelength_start.value
            < conf.supernova.luminosity_wavelength_end.value
        ):
            raise ValueError(
                "End Limit must be greater than Start Limit for Luminosity"
            )
    assert ve.type is ValueError


def test_plasma_section_config(tardis_config_verysimple):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File

    Validates:
        Initial temperature inner (must be greater than -1K)
        Initial radiative temperature (must be greater than -1K)

    Parameter
    ---------
        `tardis_config_verysimple` : YAML File

    Result
    ------
        Assertion based on validation for specified values
    """
    conf = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )
    tardis_config_verysimple["plasma"]["initial_t_inner"] = "-100 K"
    tardis_config_verysimple["plasma"]["initial_t_rad"] = "-100 K"
    with pytest.raises(ValueError) as ve:
        if (conf.plasma.initial_t_inner.value >= -1) and (
            conf.plasma.initial_t_rad.value >= -1
        ):
            raise ValueError("Initial Temperatures are Invalid")
    assert ve.type is ValueError


def test_plasma_nlte_section_config(tardis_config_verysimple_nlte):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File.

    Validates:
        nlte_ionization_species: should be included in continuum_interaction

    Parameter
    ---------
        `tardis_config_verysimple` : YAML File

    Result
    ------
        Assertion based on validation for specified values
    """
    conf = Configuration.from_config_dict(
        tardis_config_verysimple_nlte, validate=True, config_dirname="test"
    )
    tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
        "species"
    ] = [
        "He I",
    ]
    tardis_config_verysimple_nlte["plasma"]["nlte_ionization_species"] = ["H I"]
    with pytest.raises(ValueError) as ve:
        nlte_ionization_species = tardis_config_verysimple_nlte["plasma"][
            "nlte_ionization_species"
        ]

        for species in nlte_ionization_species:
            if not (
                species
                in tardis_config_verysimple_nlte["plasma"][
                    "continuum_interaction"
                ]["species"]
            ):
                raise ValueError("Nlte ionization species not in continuum.")
    assert ve.type is ValueError


def test_spectrum_section_config(tardis_config_verysimple):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File

    Validates:
        Spectrum Start & End Limits (Start < End)

    Parameter
    ---------
        `tardis_config_verysimple` : YAML File

    Result
    ------
        Assertion based on validation for specified values
    """
    conf = Configuration.from_config_dict(
        tardis_config_verysimple, validate=True, config_dirname="test"
    )
    tardis_config_verysimple["spectrum"]["start"] = "2500 angstrom"
    tardis_config_verysimple["spectrum"]["stop"] = "500 angstrom"
    with pytest.raises(ValueError) as ve:
        if not conf.spectrum.stop.value < conf.spectrum.start.value:
            raise ValueError("Start Value must be less than Stop Value")
    assert ve.type is ValueError
