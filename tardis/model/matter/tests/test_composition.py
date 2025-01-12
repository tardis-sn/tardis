import pytest
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.model.parse_geometry_configuration import parse_geometry_from_config
from tardis.io.model.parse_atom_data import parse_atom_data
from tardis.io.model.parse_composition_configuration import parse_composition_from_config


@pytest.fixture(scope="function")
def composition_instance(config_verysimple):
    time_explosion = config_verysimple.supernova.time_explosion.cgs
    geometry = parse_geometry_from_config(config_verysimple, time_explosion)
    atom_data = parse_atom_data(config_verysimple)

    composition, _ = parse_composition_from_config(
        atom_data, config_verysimple, time_explosion, geometry
    )
    return composition


def test_elemental_mass_fraction(composition_instance):
    elemental_mass_fraction = composition_instance.elemental_mass_fraction
    assert isinstance(elemental_mass_fraction, pd.DataFrame)
    np.testing.assert_allclose(elemental_mass_fraction.sum(), 1.0, rtol=1e-3)


def test_elemental_number_density(composition_instance):
    number_density = composition_instance.elemental_number_density
    assert isinstance(number_density, pd.DataFrame)
    assert np.all(number_density.values >= 0)

    # Calculating expected number density manually
    density_value = composition_instance.density.to(u.g / u.cm**3).value
    mass_fractions = composition_instance.elemental_mass_fraction
    element_masses = composition_instance.effective_element_masses.reindex(
        mass_fractions.index
    )
    
    expected_density = mass_fractions.multiply(density_value).divide(
        element_masses,
        axis=0
    )
    
    np.testing.assert_allclose(
        number_density.values,
        expected_density.values,
        rtol=1e-10
    )


@pytest.mark.parametrize("time_explosion", [10 * u.s, 100 * u.s])
def test_calculate_mass_fraction_at_time(composition_instance, time_explosion):
    if composition_instance.isotopic_mass_fraction.empty:
        result = composition_instance.calculate_mass_fraction_at_time(time_explosion)
        pd.testing.assert_frame_equal(result, composition_instance.elemental_mass_fraction)
    else:
        initial_state = composition_instance.isotopic_mass_fraction.copy()
        composition_instance.calculate_mass_fraction_at_time(time_explosion)
        
        assert not composition_instance.isotopic_mass_fraction.equals(initial_state)


def test_calculate_cell_masses(composition_instance):
    volume = 10 * u.cm**3
    density = composition_instance.density
    cell_masses = composition_instance.calculate_cell_masses(volume)

    assert isinstance(cell_masses, u.Quantity)
    assert cell_masses.unit == u.g

    expected_masses = (density * volume).to(u.g)
    np.testing.assert_allclose(cell_masses.value, expected_masses.value)


def test_calculate_elemental_cell_masses(composition_instance):
    volume = 10 * u.cm**3
    density = composition_instance.density
    elemental_masses = composition_instance.calculate_elemental_cell_masses(volume)

    assert isinstance(elemental_masses, pd.DataFrame)
    assert np.all(elemental_masses.values >= 0)

    expected_masses = composition_instance.elemental_mass_fraction * (density * volume).to(u.g).value
    np.testing.assert_allclose(elemental_masses.values, expected_masses.values)
