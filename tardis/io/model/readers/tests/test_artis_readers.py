import numpy.testing as npt
from astropy import units as u

from tardis.io.model.readers.artis import (
    read_artis_density,
    read_artis_mass_fractions,
    read_artis_model,
)


def test_artis_density_reader(example_model_file_dir):
    # Using a test ARTIS density file.
    # File: tardis_artis_density_test.dat
    time_model, velocity, mean_density = read_artis_density(
        example_model_file_dir / "artis_model.dat"
    )
    # Check that time is recognized as time
    assert time_model.unit.physical_type == "time"
    # Check velocity unit is cm/s
    assert velocity.unit == u.Unit("cm/s")
    # Dummy check for mean_density value
    npt.assert_allclose(mean_density.value[0], 10**-5.449497e-05)
    # Additional check that all density values are positive
    assert (mean_density.value > 0).all(), "All mean_density values should be positive"


def test_artis_mass_fractions_reader(example_model_file_dir):
    # Using a test ARTIS abundance file.
    mass_fractions = read_artis_mass_fractions(
        example_model_file_dir / "artis_abundances.dat"
    )
    # Verify index and data shape (adjust expected shape as appropriate)
    assert mass_fractions.index.name == "atomic_number"
    assert mass_fractions.columns.name == "cell_index"
    # Dummy check: ensure there is at least one value
    npt.assert_(mass_fractions.size > 0, "Mass fractions should not be empty")


def test_artis_model_reader(example_model_file_dir):
    # Combine both density and abundance test files.
    artis_model_data = read_artis_model(
        example_model_file_dir / "artis_model.dat",
        example_model_file_dir / "artis_abundances.dat",
    )
    # Check that the returned object has proper attributes.
    assert hasattr(artis_model_data, "time_of_model")
    assert hasattr(artis_model_data, "velocity")
    assert hasattr(artis_model_data, "mean_density")
    assert hasattr(artis_model_data, "mass_fractions")
    # Check basic units
    assert artis_model_data.time_of_model.unit.physical_type == "time"
    assert artis_model_data.velocity.unit.physical_type == "velocity"
    # Dummy check for mass fractions DataFrame (non-empty)
    assert (
        artis_model_data.mass_fractions.size > 0,
        "Mass fractions should not be empty",
    )


# ...existing code if any...
