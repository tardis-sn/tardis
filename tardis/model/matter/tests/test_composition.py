import pytest
import numpy as np
import pandas as pd
import pandas.testing as pdt
import numpy.testing as npt
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


def test_elemental_mass_fraction(test_composition_simple, regression_data):
    actual = test_composition_simple.elemental_mass_fraction
    assert isinstance(actual, pd.DataFrame)

    expected = regression_data.sync_dataframe(actual)
    pdt.assert_frame_equal(actual, expected)


def test_elemental_number_density(test_composition_simple, regression_data):
    number_density = test_composition_simple.elemental_number_density

    assert isinstance(number_density, pd.DataFrame)
    assert np.all(number_density.values >= 0)

    expected_number_density = regression_data.sync_dataframe(number_density)
    pdt.assert_frame_equal(number_density, expected_number_density)


@pytest.mark.parametrize("time_explosion", [10 * u.s, 100 * u.s])
def test_calculate_mass_fraction_at_time(test_composition_simple, time_explosion, regression_data):
    if test_composition_simple.isotopic_mass_fraction.empty:
        result = test_composition_simple.calculate_mass_fraction_at_time(time_explosion)
        expected = regression_data.sync_dataframe(result)
        pdt.assert_frame_equal(result, expected)
    else:
        initial_state = test_composition_simple.isotopic_mass_fraction.copy()
        test_composition_simple.calculate_mass_fraction_at_time(time_explosion)
        
        assert not test_composition_simple.isotopic_mass_fraction.equals(initial_state)


def test_calculate_cell_masses(test_composition_simple, regression_data):
    volume = 10 * u.cm**3
    density = test_composition_simple.density
    cell_masses = test_composition_simple.calculate_cell_masses(volume)

    expected = regression_data.sync_ndarray(cell_masses.value)
    npt.assert_allclose(cell_masses.value, expected)

    assert cell_masses.unit == u.g

    expected_masses = (density * volume).to(u.g)
    npt.assert_allclose(cell_masses.value, expected_masses.value)


def test_calculate_elemental_cell_masses(test_composition_simple, regression_data):
    volume = 10 * u.cm**3
    elemental_masses = test_composition_simple.calculate_elemental_cell_masses(volume)

    assert isinstance(elemental_masses, pd.DataFrame)
    assert np.all(elemental_masses.values >= 0)

    expected_mass_values = regression_data.sync_ndarray(elemental_masses.values)
    npt.assert_allclose(elemental_masses.values, expected_mass_values)

