from astropy import units as u
import numpy as np
import numpy.testing as nptest
import pytest
from tardis.model._generic_model import *
from tardis.io.model_reader import read_abundances_file


@pytest.fixture(scope="module")
def reference_velocity():
    """Dictionary containing reference data for `test_velocity`"""

    velocity_inner = (
        [9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500] * u.km / u.s
    ).cgs
    velocity_outer = (
        [9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000] * u.km / u.s
    ).cgs
    velocity_middle = (velocity_inner + velocity_outer) / 2
    inner_radius = velocity_inner * (2 * u.d).to("s")
    middle_radius = velocity_middle * (2 * u.d).to("s")
    outer_radius = velocity_outer * (2 * u.d).to("s")
    volume = 4 / 3 * np.pi * (outer_radius ** 3 - inner_radius ** 3)
    reference_data = {
        "velocity_inner": velocity_inner,
        "velocity_outer": velocity_outer,
        "velocity_middle": velocity_middle,
        "inner_radius": inner_radius,
        "middle_radius": middle_radius,
        "outer_radius": outer_radius,
        "volume": volume,
    }
    return reference_data


def test_property_velocity(reference_velocity):
    """Test if a correct Velocity class is created with input"""
    reference_data = reference_velocity
    velocity_generator = [i for i in range(9000, 13500, 500)]
    vel_inner, vel_outer = (
        velocity_generator[:-1] * u.km / u.s,
        velocity_generator[1:] * u.km / u.s,
    )
    velocity_class = Velocity(vel_inner, vel_outer, time=2 * u.d)

    nptest.assert_allclose(
        velocity_class.inner_radius.value, reference_data["inner_radius"].value
    )
    nptest.assert_allclose(
        velocity_class.middle_radius.value,
        reference_data["middle_radius"].value,
    )
    nptest.assert_allclose(
        velocity_class.outer_radius.value, reference_data["outer_radius"].value
    )
    nptest.assert_allclose(
        velocity_class.inner.value, reference_data["velocity_inner"].value
    )
    nptest.assert_allclose(
        velocity_class.middle.value, reference_data["velocity_middle"].value
    )
    nptest.assert_allclose(
        velocity_class.outer.value, reference_data["velocity_outer"].value
    )
    nptest.assert_allclose(
        velocity_class.volume.value,
        reference_data["volume"].value,
        rtol=1.0e-15,
    )
    assert velocity_class.number_of_cells == 8


@pytest.fixture(scope="module")
def reference_density():
    """Reference density for `test_property_density"""
    initial_mass_density = (
        [
            6.2419e-10,
            6.1502e-10,
            6.0492e-10,
            5.9475e-10,
            5.8468e-10,
            5.7473e-10,
            5.6493e-10,
            5.5529e-10,
        ]
        * u.g
        / u.cm ** 3
    )

    density_after_time = (
        initial_mass_density * ((2 * u.d).to("s") / 1 * u.s) ** -3
    )
    #    electron_density =
    reference_data = {
        "initial_mass_density": initial_mass_density,
        "density_after_time": density_after_time,
    }

    return reference_data


def test_property_density(reference_density):
    """Test if a correct Density class is created with input"""
    reference_data = reference_density
    density = (
        [
            6.2419e-10,
            6.1502e-10,
            6.0492e-10,
            5.9475e-10,
            5.8468e-10,
            5.7473e-10,
            5.6493e-10,
            5.5529e-10,
        ]
        * u.g
        / u.cm ** 3
    )
    density_class = Density(density, time=2 * u.d)
    nptest.assert_allclose(
        density_class.mass.value, reference_data["density_after_time"].value
    )
    nptest.assert_allclose(
        density_class.mass_0.value, reference_data["initial_mass_density"].value
    )
    assert density_class.number_of_cells == 8


@pytest.fixture(scope="module")
def reference_abundance():
    """Read reference abundance for `test_property_abundances`"""
    abundance = np.array(
        [
            [0.0, 0.4, 0.35, 0.35, 0.38, 0.4, 0.38, 0.1],
            [0.99, 0.58, 0.6, 0.6, 0.57, 0.58, 0.57, 0.8],
        ]
    )
    isotope_abundance = np.array(
        [[0.01, 0.02, 0.05, 0.05, 0.05, 0.02, 0.05, 0.1]]
    )
    decayed = np.array(
        [
            [0.0, 0.4, 0.35, 0.35, 0.38, 0.4, 0.38, 0.1],
            [0.99, 0.58, 0.6, 0.6, 0.57, 0.58, 0.57, 0.8],
            [
                1.889189267968716e-05,
                3.778378535937432e-05,
                9.44594633984358e-05,
                9.44594633984358e-05,
                9.44594633984358e-05,
                3.778378535937432e-05,
                9.44594633984358e-05,
                0.0001889189267968716,
            ],
            [
                0.0020214305565343163,
                0.004042861113068633,
                0.010107152782671582,
                0.010107152782671582,
                0.010107152782671582,
                0.004042861113068633,
                0.010107152782671582,
                0.020214305565343163,
            ],
            [
                0.007959677550785997,
                0.015919355101571993,
                0.039798387753929985,
                0.039798387753929985,
                0.039798387753929985,
                0.015919355101571993,
                0.039798387753929985,
                0.07959677550785997,
            ],
        ]
    )
    reference_data = {
        "abundance": abundance,
        "isotope_abundance": isotope_abundance,
        "decayed_abundance": decayed,
    }
    return reference_data


def test_property_abundances(reference_abundance):
    reference_data = reference_abundance
    abundance_filename = "data/custom_abundance_tardis_gm.dat"
    index, abundance, isotope_abundance = read_abundances_file(
        abundance_filename, "custom_composition"
    )
    abundance_class = Abundances(abundance, isotope_abundance, time=2 * u.d)
    nptest.assert_allclose(
        abundance_class.elemental_0.to_numpy(), reference_data["abundance"]
    )
    nptest.assert_allclose(
        abundance_class.isotope_0.to_numpy(),
        reference_data["isotope_abundance"],
    )
    nptest.assert_allclose(
        abundance_class.elemental.T.to_numpy(),
        reference_data["decayed_abundance"],
        rtol=1e-3,
    )


@pytest.fixture(scope="module")
def reference_radiation_field():
    radiative_temperature = ([7000] * 8) * u.K
    dilution_factor = [0.3] * 8
    reference_data = {
        "radiative_temperature": radiative_temperature,
        "dilution_factor": dilution_factor,
    }
    return reference_data


def test_property_radiation_field(reference_radiation_field):
    reference_data = reference_radiation_field
    radiative_temperature = [
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
    ] * u.K
    dilution_factor = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    radiation_field = RadiationField(
        radiative_temperature, dilution_factor, time=2 * u.d
    )
    nptest.assert_allclose(
        radiation_field.radiative_temperature,
        reference_data["radiative_temperature"],
    )
    nptest.assert_allclose(
        radiation_field.dilution_factor, reference_data["dilution_factor"]
    )


def test_generic_model(
    reference_velocity,
    reference_density,
    reference_abundance,
    reference_radiation_field,
):
    reference_velocity = reference_velocity
    reference_density = reference_density
    reference_abundance = reference_abundance
    reference_radiation_field = reference_radiation_field
    velocity_generator = [i for i in range(9000, 13500, 500)]
    vel_inner, vel_outer = (
        velocity_generator[:-1] * u.km / u.s,
        velocity_generator[1:] * u.km / u.s,
    )
    density = (
        [
            6.2419e-10,
            6.1502e-10,
            6.0492e-10,
            5.9475e-10,
            5.8468e-10,
            5.7473e-10,
            5.6493e-10,
            5.5529e-10,
        ]
        * u.g
        / u.cm ** 3
    )
    abundance_filename = "data/custom_abundance_tardis_gm.dat"
    index, abundance, isotope_abundance = read_abundances_file(
        abundance_filename, "custom_composition"
    )
    radiative_temperature = [
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
        7000,
    ] * u.K
    dilution_factor = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]

    velocity_class = Velocity(vel_inner, vel_outer, time=2 * u.d)
    density_class = Density(density, time=1 * u.d)
    abundance_class = Abundances(abundance, isotope_abundance, time=1.7 * u.d)
    radiation_field = RadiationField(
        radiative_temperature, dilution_factor, time=1 * u.d
    )

    generic_model = GenericModel(
        velocity_class, density_class, abundance_class, radiation_field
    )

    for prop in generic_model:
        assert getattr(generic_model, prop).time == 2 * u.d
        assert getattr(generic_model, prop).number_of_cells == 8

    nptest.assert_allclose(
        generic_model.velocity.inner_radius.value,
        reference_velocity["inner_radius"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.middle_radius.value,
        reference_velocity["middle_radius"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.outer_radius.value,
        reference_velocity["outer_radius"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.inner.value,
        reference_velocity["velocity_inner"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.middle.value,
        reference_velocity["velocity_middle"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.outer.value,
        reference_velocity["velocity_outer"].value,
    )
    nptest.assert_allclose(
        generic_model.velocity.volume.value,
        reference_velocity["volume"].value,
        rtol=1.0e-15,
    )
    nptest.assert_allclose(
        generic_model.density.mass.value,
        reference_density["density_after_time"].value,
    )
    nptest.assert_allclose(
        generic_model.density.mass_0.value,
        reference_density["initial_mass_density"].value,
    )
    nptest.assert_allclose(
        generic_model.abundance.elemental_0.to_numpy(),
        reference_abundance["abundance"],
    )
    nptest.assert_allclose(
        generic_model.abundance.isotope_0.to_numpy(),
        reference_abundance["isotope_abundance"],
    )
    nptest.assert_allclose(
        generic_model.abundance.elemental.T.to_numpy(),
        reference_abundance["decayed_abundance"],
        rtol=1e-3,
    )
    nptest.assert_allclose(
        generic_model.radiation_field.radiative_temperature,
        reference_radiation_field["radiative_temperature"],
    )
    nptest.assert_allclose(
        generic_model.radiation_field.dilution_factor,
        reference_radiation_field["dilution_factor"],
    )
