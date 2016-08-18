import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy import constants as const
from tardis.io.atomic import atomic_number2symbol, symbol2atomic_number


@pytest.fixture
def basic_atom_data_he(atom_data_he):
    return atom_data_he.atom_data


@pytest.fixture
def ionization_data_he(atom_data_he):
    return atom_data_he.ionization_data


@pytest.fixture
def levels_he(atom_data_he):
    return atom_data_he.levels


@pytest.fixture
def lines_he(atom_data_he):
    return atom_data_he.lines


def test_atom_data_basic_atom_data(basic_atom_data_he):
    assert basic_atom_data_he.loc[2, 'symbol'] == 'He'
    assert_quantity_allclose(
        basic_atom_data_he.loc[2, 'mass'] * u.Unit('g'),
        4.002602 * const.u.cgs
     )


def test_atom_data_ionization_data(ionization_data_he):
    assert_quantity_allclose(
        ionization_data_he.loc[(2, 1), 'ionization_energy'] * u.Unit('erg'),
        24.587387936 * u.Unit('eV')
    )


def test_atom_data_levels(levels_he):
    levels_he = levels_he.set_index(['atomic_number', 'ion_number', 'level_number'])
    assert_quantity_allclose(
        u.Quantity(levels_he.loc[(2, 0, 2), 'energy'], u.Unit('erg')).to(u.Unit('cm-1'), equivalencies=u.spectral()),
        166277.542 * u.Unit('cm-1')
    )


def test_atom_data_lines(lines_he):
    lines_he = lines_he.set_index(['atomic_number', 'ion_number',
                                   'level_number_lower', 'level_number_upper'])

    assert_quantity_allclose(
        lines_he.loc[(2,0,0,6), 'wavelength_cm'] * u.Unit('cm'),
        584.335 * u.Unit('Angstrom')
    )


def test_atomic_symbol():
    assert atomic_number2symbol[14] == 'Si'


def test_atomic_symbol_reverse():
    assert symbol2atomic_number['Si'] == 14
