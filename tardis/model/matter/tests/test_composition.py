import pytest
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.model.parse_geometry_configuration import parse_geometry_from_config
from tardis.io.model.parse_atom_data import parse_atom_data
from tardis.io.model.parse_composition_configuration import parse_composition_from_config


@pytest.fixture(scope="module")
def test_composition_simple(config_verysimple, atomic_dataset):
    time_explosion = config_verysimple.supernova.time_explosion.cgs
    geometry = parse_geometry_from_config(config_verysimple, time_explosion)

    composition, _ = parse_composition_from_config(
        atomic_dataset, config_verysimple, time_explosion, geometry
    )
    return composition


def test_elemental_mass_fraction(test_composition_simple):
    elemental_mass_fraction = test_composition_simple.elemental_mass_fraction
    assert isinstance(elemental_mass_fraction, pd.DataFrame)
    np.testing.assert_allclose(elemental_mass_fraction.sum(), 1.0)


def test_elemental_number_density(test_composition_simple):
    number_density = test_composition_simple.elemental_number_density
    assert isinstance(number_density, pd.DataFrame)
    assert np.all(number_density.values >= 0)

    # Calculating expected number density manually
    density_value = test_composition_simple.density.to(u.g / u.cm**3).value
    mass_fractions = test_composition_simple.elemental_mass_fraction
    element_masses = test_composition_simple.effective_element_masses.reindex(
        mass_fractions.index
    )
    
    expected_number_density = mass_fractions.multiply(density_value).divide(
        element_masses,
        axis=0
    )
    
    np.testing.assert_allclose(
        number_density.values,
        expected_number_density.values,
        rtol=1e-10
    )


@pytest.mark.parametrize("time_explosion", [10 * u.s, 100 * u.s])
def test_calculate_mass_fraction_at_time(test_composition_simple, time_explosion):
    if test_composition_simple.isotopic_mass_fraction.empty:
        result = test_composition_simple.calculate_mass_fraction_at_time(time_explosion)
        pd.testing.assert_frame_equal(result, test_composition_simple.elemental_mass_fraction)
    else:
        initial_state = test_composition_simple.isotopic_mass_fraction.copy()
        test_composition_simple.calculate_mass_fraction_at_time(time_explosion)
        
        assert not test_composition_simple.isotopic_mass_fraction.equals(initial_state)


def test_calculate_cell_masses(test_composition_simple):
    volume = 10 * u.cm**3
    density = test_composition_simple.density
    cell_masses = test_composition_simple.calculate_cell_masses(volume)

    assert isinstance(cell_masses, u.Quantity)
    assert cell_masses.unit == u.g

    expected_masses = (density * volume).to(u.g)
    np.testing.assert_allclose(cell_masses.value, expected_masses.value)


def test_calculate_elemental_cell_masses(test_composition_simple):
    volume = 10 * u.cm**3
    density = test_composition_simple.density
    elemental_masses = test_composition_simple.calculate_elemental_cell_masses(volume)

    assert isinstance(elemental_masses, pd.DataFrame)
    assert np.all(elemental_masses.values >= 0)

    expected_masses = test_composition_simple.elemental_mass_fraction * (density * volume).to(u.g).value
    np.testing.assert_allclose(elemental_masses.values, expected_masses.values)

