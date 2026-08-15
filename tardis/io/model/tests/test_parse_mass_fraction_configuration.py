from pathlib import Path

import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.parse_geometry_configuration import (
    parse_geometry_from_config,
)
from tardis.io.model.parse_mass_fraction_configuration import (
    parse_mass_fractions_from_config,
)
from tardis.model.geometry.radial1d_homologous import (
    HomologousRadial1DGeometry,
)

HYDROGEN = (1, -1)
HELIUM = (2, -1)
NI56 = (28, 56)
CO56 = (27, 56)
FE56 = (26, 56)


@pytest.fixture
def config_dict_abundances(example_configuration_dir: Path) -> Configuration:
    """Config whose abundance section is replaced by a plain dict, the way a
    model is built in a notebook without going through the schema validator.
    """
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    config.model.abundances = {
        "type": "uniform",
        "H": 0.70,
        "He": 0.25,
        "Ni56": 0.05,
    }
    return config


@pytest.fixture
def geometry_dict_abundances(
    config_dict_abundances: Configuration,
) -> HomologousRadial1DGeometry:
    return parse_geometry_from_config(
        config_dict_abundances,
        config_dict_abundances.supernova.time_explosion.cgs,
    )


def test_dict_abundances_without_model_isotope_time_0(
    config_dict_abundances: Configuration,
    geometry_dict_abundances: HomologousRadial1DGeometry,
) -> None:
    """A dict abundance section without ``model_isotope_time_0`` falls back to
    the schema default of nan seconds, which leaves the isotopes undecayed.
    """
    nuclide_mass_fractions = parse_mass_fractions_from_config(
        config_dict_abundances,
        geometry_dict_abundances,
        config_dict_abundances.supernova.time_explosion.cgs,
    )

    no_of_shells = geometry_dict_abundances.no_of_shells
    assert nuclide_mass_fractions.shape[1] == no_of_shells

    npt.assert_allclose(
        nuclide_mass_fractions.loc[HYDROGEN, :].to_numpy(),
        [0.70] * no_of_shells,
    )
    npt.assert_allclose(
        nuclide_mass_fractions.loc[HELIUM, :].to_numpy(),
        [0.25] * no_of_shells,
    )
    npt.assert_allclose(
        nuclide_mass_fractions.loc[NI56, :].to_numpy(),
        [0.05] * no_of_shells,
    )
    assert CO56 not in nuclide_mass_fractions.index


def test_dict_abundances_with_model_isotope_time_0(
    config_dict_abundances: Configuration,
    geometry_dict_abundances: HomologousRadial1DGeometry,
) -> None:
    """An explicit ``model_isotope_time_0`` still decays the isotopes from
    that time to the time of the explosion.
    """
    config_dict_abundances.model.abundances["model_isotope_time_0"] = 0 * u.s

    nuclide_mass_fractions = parse_mass_fractions_from_config(
        config_dict_abundances,
        geometry_dict_abundances,
        config_dict_abundances.supernova.time_explosion.cgs,
    )

    no_of_shells = geometry_dict_abundances.no_of_shells

    npt.assert_allclose(
        nuclide_mass_fractions.loc[HYDROGEN, :].to_numpy(),
        [0.70] * no_of_shells,
    )
    assert (nuclide_mass_fractions.loc[NI56, :] < 0.05).all()
    assert (nuclide_mass_fractions.loc[CO56, :] > 0).all()
    assert FE56 in nuclide_mass_fractions.index

    # the mass 56 decay chain conserves mass
    chain_56 = nuclide_mass_fractions.loc[[NI56, CO56, FE56], :].sum(axis=0)
    npt.assert_allclose(chain_56.to_numpy(), [0.05] * no_of_shells)


def test_validated_config_abundances_are_unchanged(
    example_configuration_dir: Path,
) -> None:
    """Reading the key through the dict interface returns the schema value for
    a config that did go through the validator.
    """
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    geometry = parse_geometry_from_config(
        config, config.supernova.time_explosion.cgs
    )

    nuclide_mass_fractions = parse_mass_fractions_from_config(
        config, geometry, config.supernova.time_explosion.cgs
    )

    expected_mass_fractions = {
        (8, -1): 0.19,
        (12, -1): 0.03,
        (14, -1): 0.52,
        (16, -1): 0.19,
        (18, -1): 0.04,
        (20, -1): 0.03,
    }
    for nuclide_index, mass_fraction in expected_mass_fractions.items():
        npt.assert_allclose(
            nuclide_mass_fractions.loc[nuclide_index, :].to_numpy(),
            [mass_fraction] * geometry.no_of_shells,
        )
