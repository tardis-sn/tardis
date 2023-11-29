# tests for the config reader module
import os
from attr import validate
import pytest
import pandas as pd
from numpy.testing import assert_almost_equal
from jsonschema.exceptions import ValidationError

from tardis.io.configuration import config_reader
from astropy.units import Quantity
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.standard_plasmas import assemble_plasma


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

    parsed_convergence_section = (
        config_reader.Configuration.parse_convergence_section(
            test_convergence_section
        )
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
    ] = Quantity("2.0e4 km/s")
    tardis_config_verysimple["model"]["structure"]["velocity"][
        "stop"
    ] = Quantity("1.1e4 km/s")

    with pytest.raises(ValueError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )


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
    tardis_config_verysimple["supernova"]["time_explosion"] = Quantity(
        "-10 day"
    )
    with pytest.raises(ValueError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )

    tardis_config_verysimple["supernova"]["time_explosion"] = Quantity("10 day")
    tardis_config_verysimple["supernova"][
        "luminosity_wavelength_start"
    ] = Quantity("15 angstrom")
    tardis_config_verysimple["supernova"][
        "luminosity_wavelength_end"
    ] = Quantity("0 angstrom")
    with pytest.raises(ValueError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )


@pytest.mark.parametrize("key", ["initial_t_inner", "initial_t_rad"])
def test_plasma_section_config(key, tardis_config_verysimple):
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
    tardis_config_verysimple["plasma"][key] = Quantity("-100 K")
    with pytest.raises(ValueError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )


def test_plasma_nlte_section_root_config(
    tardis_config_verysimple_nlte,
    nlte_raw_model_root,
    nlte_atom_data,
):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File.

    Validates:
        nlte_ionization_species: should be included in continuum_interaction

    Parameter
    ---------
        `tardis_config_verysimple_nlte_root` : YAML File
        `nlte_raw_model` : A simple model
        `nlte_atom_data` : An example atomic dataset

    Result
    ------
        Assertion based on validation for specified values
    """
    tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
        "species"
    ] = [
        "He I",
    ]
    tardis_config_verysimple_nlte["plasma"]["nlte_ionization_species"] = ["H I"]
    config = Configuration.from_config_dict(tardis_config_verysimple_nlte)
    with pytest.raises(PlasmaConfigError) as ve:
        assemble_plasma(config, nlte_raw_model_root, nlte_atom_data)


def test_plasma_nlte_section_lu_config(
    tardis_config_verysimple_nlte,
    nlte_raw_model_lu,
    nlte_atom_data,
):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File.

    Validates:
        nlte_ionization_species: should be included in continuum_interaction

    Parameter
    ---------
        `tardis_config_verysimple_nlte_root` : YAML File
        `nlte_raw_model` : A simple model
        `nlte_atom_data` : An example atomic dataset

    Result
    ------
        Assertion based on validation for specified values
    """
    tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
        "species"
    ] = [
        "He I",
    ]
    tardis_config_verysimple_nlte["plasma"]["nlte_ionization_species"] = ["H I"]
    tardis_config_verysimple_nlte["plasma"]["nlte_solver"] = "lu"
    config = Configuration.from_config_dict(tardis_config_verysimple_nlte)
    with pytest.raises(PlasmaConfigError) as ve:
        assemble_plasma(config, nlte_raw_model_lu, nlte_atom_data)


def test_plasma_nlte_root_exc_section_config(
    tardis_config_verysimple_nlte, nlte_raw_model_root, nlte_atom_data
):
    """
    Configuration Validation Test for Plasma Section of the Tardis Config YAML File.

    Validates:
        nlte_excitation_species: should be included in continuum_interaction

    Parameter
    ---------
        `tardis_config_verysimple_nlte_root` : YAML File
        `nlte_raw_model` : A simple model
        `nlte_atom_data` : An example atomic dataset

    Result
    ------
        Assertion based on validation for specified values
    """
    tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
        "species"
    ] = [
        "He I",
    ]
    tardis_config_verysimple_nlte["plasma"]["nlte_excitation_species"] = ["H I"]
    tardis_config_verysimple_nlte["plasma"]["nlte_solver"] = "root"
    config = Configuration.from_config_dict(tardis_config_verysimple_nlte)
    with pytest.raises(PlasmaConfigError):
        plasma = assemble_plasma(config, nlte_raw_model_root, nlte_atom_data)


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
    tardis_config_verysimple["spectrum"]["start"] = Quantity("2500 angstrom")
    tardis_config_verysimple["spectrum"]["stop"] = Quantity("500 angstrom")
    with pytest.raises(ValueError):
        conf = Configuration.from_config_dict(
            tardis_config_verysimple, validate=True, config_dirname="test"
        )
