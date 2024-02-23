import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from tardis import constants as const


@pytest.fixture
def basic_atom_data(kurucz_atomic_data):
    return kurucz_atomic_data.atom_data


@pytest.fixture
def ionization_data(kurucz_atomic_data):
    return kurucz_atomic_data.ionization_data


@pytest.fixture
def levels(kurucz_atomic_data):
    return kurucz_atomic_data.levels


@pytest.fixture
def lines(kurucz_atomic_data):
    return kurucz_atomic_data.lines


def test_atom_data_basic_atom_data(basic_atom_data):
    assert basic_atom_data.loc[2, "symbol"] == "He"
    assert_quantity_allclose(
        basic_atom_data.at[2, "mass"] * u.Unit("g"), 4.002602 * const.u.cgs
    )


def test_atom_data_ionization_data(ionization_data):
    assert_quantity_allclose(
        ionization_data.loc[(2, 1)] * u.Unit("erg"), 24.587387936 * u.Unit("eV")
    )


def test_atom_data_levels(levels):
    assert_quantity_allclose(
        u.Quantity(levels.at[(2, 0, 2), "energy"], u.Unit("erg")).to(
            u.Unit("cm-1"), equivalencies=u.spectral()
        ),
        166277.542 * u.Unit("cm-1"),
    )


def test_atom_data_lines(lines):
    assert_quantity_allclose(
        lines.loc[(2, 0, 0, 6), "wavelength_cm"].values[0] * u.Unit("cm"),
        584.335 * u.Unit("Angstrom"),
    )


def test_atomic_reprepare(kurucz_atomic_data):
    kurucz_atomic_data.prepare_atom_data(
        [14, 20],
        line_interaction_type="scatter",
        nlte_species=[],
        continuum_interaction_species=[],
    )
    lines = kurucz_atomic_data.lines.reset_index()
    assert lines["atomic_number"].isin([14, 20]).all()
    assert len(lines.loc[lines["atomic_number"] == 14]) > 0
    assert len(lines.loc[lines["atomic_number"] == 20]) > 0
