#tests for util module

import pytest
from astropy import units as u
import numpy as np
from tardis import atomic
from tardis.util import species_string_to_tuple, species_tuple_to_string

from tardis.util import MalformedSpeciesError, MalformedElementSymbolError, MalformedQuantityError
from tardis.util import (int_to_roman, parse_quantity,
                         element_symbol2atomic_number, atomic_number2element_symbol,
                         reformat_element_symbol, quantity_linspace)


def test_malformed_species_error():
    malformed_species_error = MalformedSpeciesError('He')
    assert malformed_species_error.malformed_element_symbol == 'He'
    assert str(malformed_species_error) == 'Expecting a species notation (e.g. "Si 2", "Si II", "Fe IV") - supplied He'


def test_malformed_elements_symbol_error():
    malformed_elements_symbol_error = MalformedElementSymbolError('Hx')
    assert malformed_elements_symbol_error.malformed_element_symbol == 'Hx'
    assert str(malformed_elements_symbol_error) == 'Expecting an atomic symbol (e.g. Fe) - supplied Hx'


def test_malformed_quantity_error():
    malformed_quantity_error = MalformedQuantityError('abcd')
    assert malformed_quantity_error.malformed_quantity_string == 'abcd'
    assert str(malformed_quantity_error) == 'Expecting a quantity string(e.g. "5 km/s") for keyword - supplied abcd'


@pytest.mark.parametrize(['test_input', 'expected_result'], [(1, 'I'), (5, 'V'), (19, 'XIX'), (556, 'DLVI'),
                                                             (1400, 'MCD'), (1999, 'MCMXCIX'), (3000, 'MMM')])
def test_int_to_roman(test_input, expected_result):
    assert int_to_roman(test_input) == expected_result

    with pytest.raises(TypeError ): int_to_roman(1.5)
    with pytest.raises(ValueError): int_to_roman(0)
    with pytest.raises(ValueError): int_to_roman(4000)


def test_parse_quantity():
    q1 = parse_quantity('5 km/s')
    assert q1.value == 5.
    assert q1.unit == u.Unit('km/s')

    with pytest.raises(MalformedQuantityError):
        parse_quantity(5)

    with pytest.raises(MalformedQuantityError):
        parse_quantity('abcd')

    with pytest.raises(MalformedQuantityError):
        parse_quantity('a abcd')

    with pytest.raises(MalformedQuantityError):
        parse_quantity('5 abcd')


def test_atomic_number2element_symbol():
    assert atomic_number2element_symbol(14) == 'Si'


@pytest.mark.parametrize("element_symbol, atomic_number", [
    ('sI', 14),
    ('ca', 20),
    ('Fe', 26)
])
def test_element_symbol2atomic_number(element_symbol, atomic_number):
    assert element_symbol2atomic_number(element_symbol) == atomic_number

    with pytest.raises(MalformedElementSymbolError):
        element_symbol2atomic_number('Hx')


@pytest.mark.parametrize("unformatted_element_string, formatted_element_string", [
    ('si', 'Si'),
    ('sI', 'Si'),
    ('Si', 'Si'),
    ('c', 'C'),
    ('C', 'C'),
])
def test_element_symbol_reformatter(unformatted_element_string, formatted_element_string):
    assert reformat_element_symbol(unformatted_element_string) == formatted_element_string


@pytest.mark.parametrize("species_string, species_tuple", [
    ('si ii', (14, 1)),
    ('si 2', (14, 1)),
    ('si ix', (14, 8)),
])
def test_species_string_to_tuple(species_string, species_tuple):
    assert species_string_to_tuple(species_string) == species_tuple

    with pytest.raises(MalformedSpeciesError):
        species_string_to_tuple('II')

    with pytest.raises(MalformedSpeciesError):
        species_string_to_tuple('He Si')

    with pytest.raises(ValueError):
        species_string_to_tuple('He IX')


@pytest.mark.parametrize(['species_tuple', 'roman_numerals', 'species_string'], [
    ((14, 1), True, 'Si II'), ((14, 3), True, 'Si IV'), ((14, 8), True, 'Si IX'),
    ((14, 1), False, 'Si 1'), ((14, 3), False, 'Si 3'), ((14, 8), False, 'Si 8'),
])
def test_species_tuple_to_string(species_tuple, roman_numerals, species_string):
    assert species_tuple_to_string(species_tuple, roman_numerals=roman_numerals) == species_string


@pytest.mark.parametrize(['start', 'stop', 'num', 'expected'], [
    (u.Quantity(1, 'km/s'), u.Quantity(5, 'km/s'), 5, u.Quantity(np.array([1., 2., 3., 4., 5.]), 'km/s')),
    (u.Quantity(0.5, 'eV'), u.Quantity(0.6, 'eV'), 3, u.Quantity(np.array([0.5, 0.55, 0.6]), 'eV'))
])
def test_quantity_linspace(start, stop, num, expected):
    obtained = quantity_linspace(start, stop, num)
    assert obtained.unit == expected.unit
    assert obtained.value.all() == expected.value.all()

    with pytest.raises(ValueError):
        quantity_linspace(u.Quantity(0.5, 'eV'), '0.6 eV', 3)
