from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.artis.readers import (
    read_artis_composition,
    read_artis_density,
    read_artis_mass_fractions,
    read_artis_model,
)
from tardis.io.model.parse_composition_configuration import (
    parse_composition_from_config,
)
from tardis.io.model.parse_geometry_configuration import (
    parse_geometry_from_config,
)
from tardis.io.model.readers.base import (
    read_density_file,
    read_mass_fractions_file,
)


@pytest.fixture
def artis_data_dir():
    return Path("tardis/io/model/artis/tests/data")


@pytest.fixture
def artis_density_fname(artis_data_dir):
    return artis_data_dir / "artis_model.dat"


def test_artis_density_reader(artis_density_fname: str):
    # Using a test ARTIS density file.
    # File: tardis_artis_density_test.dat
    time_model, velocity, mean_density, isotopic_mass_fractions = (
        read_artis_density(artis_density_fname, legacy_return=False)
    )
    # Check that time is recognized as time
    assert time_model.unit.physical_type == "time"
    # Check velocity unit is cm/s
    assert velocity.unit == u.Unit("cm/s")
    assert len(mean_density) == len(velocity) - 1
    assert isotopic_mass_fractions.shape[1] == len(mean_density)
    assert isotopic_mass_fractions.columns[0] == 2
    # The first ARTIS row is the unused inner shell.
    npt.assert_allclose(mean_density.value[0], 10**3.013383e-04)
    # Additional check that all density values are positive
    assert (mean_density.value > 0).all(), (
        "All mean_density values should be positive"
    )


def test_artis_mass_fractions_reader(artis_data_dir):
    # Using a test ARTIS abundance file.
    mass_fractions = read_artis_mass_fractions(
        artis_data_dir / "artis_abundances.dat"
    )
    # Verify index and data shape (adjust expected shape as appropriate)
    assert mass_fractions.index.name == "atomic_number"
    assert mass_fractions.columns.name == "cell_index"
    assert mass_fractions.shape == (30, 69)
    assert mass_fractions.columns[0] == 2
    assert mass_fractions.columns[-1] == 70
    # Dummy check: ensure there is at least one value
    assert mass_fractions.size > 0, "Mass fractions should not be empty"


def test_artis_model_reader(artis_data_dir):
    # Combine both density and abundance test files.
    artis_model_data = read_artis_model(
        artis_data_dir / "artis_model.dat",
        artis_data_dir / "artis_abundances.dat",
    )
    # Check that the returned object has proper attributes.
    assert hasattr(artis_model_data, "time_of_model")
    assert hasattr(artis_model_data, "velocity")
    assert hasattr(artis_model_data, "mean_density")
    assert hasattr(artis_model_data, "mass_fractions")
    # Check basic units
    assert artis_model_data.time_of_model.unit.physical_type == "time"
    assert artis_model_data.velocity.unit.physical_type == "velocity"
    assert (
        len(artis_model_data.mean_density) == len(artis_model_data.velocity) - 1
    )
    assert artis_model_data.mass_fractions.shape[1] == len(
        artis_model_data.mean_density
    )
    assert artis_model_data.isotope_mass_fractions.shape == (4, 69)
    npt.assert_allclose(
        artis_model_data.isotope_mass_fractions.loc[(28, 56), 2],
        0.9788986,
    )
    npt.assert_allclose(artis_model_data.mass_fractions.sum(axis=0), 1.0)
    # Dummy check for mass fractions DataFrame (non-empty)
    assert artis_model_data.mass_fractions.size > 0, (
        "Mass fractions should not be empty"
    )


def test_simple_legacy_read_artis_density(artis_density_fname: str):
    time_of_model, velocity, mean_density = read_artis_density(
        artis_density_fname, legacy_return=True
    )

    assert np.isclose(0.00114661 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(
        mean_density[23],
        0.2250048 * u.g / u.cm**3,
        atol=1.0e-6 * u.g / u.cm**3,
    )
    assert len(mean_density) == 69
    assert len(velocity) == len(mean_density) + 1


def test_artis_generic_readers_shell_slicing(artis_data_dir):
    time_of_model, velocity, mean_density, electron_densities, temperature = (
        read_density_file(artis_data_dir / "artis_model.dat", "artis")
    )
    index, mass_fractions, isotope_mass_fractions = read_mass_fractions_file(
        artis_data_dir / "artis_abundances.dat",
        "artis",
        density_filename=artis_data_dir / "artis_model.dat",
    )

    assert time_of_model.unit.physical_type == "time"
    assert electron_densities is None
    assert temperature is None
    assert len(mean_density) == len(velocity) - 1
    assert len(index) == len(mean_density)
    assert mass_fractions.shape[1] == len(mean_density)
    assert isotope_mass_fractions.shape == (4, len(mean_density))
    assert list(isotope_mass_fractions.index) == [
        (28, 56),
        (27, 56),
        (26, 52),
        (24, 48),
    ]
    direct_elemental_mass_fractions = read_artis_mass_fractions(
        artis_data_dir / "artis_abundances.dat"
    )
    reconstructed_elemental_mass_fractions = mass_fractions.add(
        isotope_mass_fractions.groupby(level="atomic_number").sum(),
        fill_value=0.0,
    )
    npt.assert_allclose(
        reconstructed_elemental_mass_fractions.to_numpy(),
        direct_elemental_mass_fractions.to_numpy(),
    )
    npt.assert_allclose(
        isotope_mass_fractions.loc[(28, 56), 0],
        0.9788986,
    )


def test_artis_composition_preserves_explicit_isotopes(artis_data_dir):
    _, elemental_mass_fractions, isotope_mass_fractions = (
        read_artis_composition(
            artis_data_dir / "artis_model.dat",
            artis_data_dir / "artis_abundances.dat",
        )
    )
    original_elemental_mass_fractions = read_artis_mass_fractions(
        artis_data_dir / "artis_abundances.dat"
    )

    npt.assert_allclose(
        isotope_mass_fractions.loc[(28, 56), 2],
        0.9788986,
    )
    npt.assert_allclose(
        elemental_mass_fractions.loc[28, 2]
        + isotope_mass_fractions.loc[(28, 56), 2],
        original_elemental_mass_fractions.loc[28, 2],
    )
    npt.assert_allclose(
        elemental_mass_fractions.sum(axis=0)
        + isotope_mass_fractions.sum(axis=0),
        1.0,
    )


def test_artis_isotopes_reach_composition(artis_data_dir):
    config = Configuration.from_yaml(
        artis_data_dir / "tardis_configv1_artis_density.yml"
    )
    config.model.abundances.type = "file"
    config.model.abundances.filename = "artis_abundances.dat"
    config.model.abundances.filetype = "artis"
    time_explosion = config.supernova.time_explosion.cgs
    geometry = parse_geometry_from_config(config, time_explosion)

    composition, electron_densities = parse_composition_from_config(
        None, config, time_explosion, geometry
    )

    expected_isotopes = [(28, 56), (27, 56), (26, 52), (24, 48)]
    assert electron_densities is None
    assert list(composition.isotopic_mass_fraction.index) == expected_isotopes
    assert list(composition.isotope_masses.index) == expected_isotopes
    assert composition.isotopic_mass_fraction.shape[1] == len(
        composition.density
    )
    assert composition.isotope_masses.shape == (
        composition.isotopic_mass_fraction.shape
    )
    assert list(composition.isotope_masses.columns) == list(
        composition.isotopic_mass_fraction.columns
    )
    assert (composition.nuclide_mass_fraction.values >= 0).all()
    npt.assert_allclose(
        composition.isotopic_mass_fraction.loc[(28, 56), 0],
        0.9788986,
    )
