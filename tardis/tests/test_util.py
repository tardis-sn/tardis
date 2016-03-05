#tests for util module

import pytest
from astropy import units as u
from tardis import atomic
from tardis.util import species_string_to_tuple, species_tuple_to_string, parse_quantity, element_symbol2atomic_number, atomic_number2element_symbol, reformat_element_symbol, MalformedQuantityError

from tardis.util import (MalformedSpeciesError, MalformedElementSymbolError)


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


def test_quantity_parser_normal():
    q1 = parse_quantity('5 km/s')
    assert q1.value == 5.
    assert q1.unit == u.Unit('km/s')

def test_quantity_parser_malformed_quantity1():
    with pytest.raises(MalformedQuantityError):
        q1 = parse_quantity('abcd')

def test_quantity_parser_malformed_quantity2():
    with pytest.raises(MalformedQuantityError):
        q1 = parse_quantity('5 abcd')

def test_atomic_number2element_symbol():
    assert atomic_number2element_symbol(14) == 'Si'


@pytest.mark.parametrize("element_symbol, atomic_number", [
    ('sI', 14),
    ('ca', 20),
    ('Fe', 26)
])
def test_element_symbol2atomic_number(element_symbol, atomic_number):
    assert element_symbol2atomic_number(element_symbol) == atomic_number


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
def test_species_string_to_species_tuple(species_string, species_tuple):
    assert species_string_to_tuple(species_string) == species_tuple


@pytest.mark.parametrize("species_string, species_tuple", [
    ('Si II', (14, 1)),
    ('Si IV', (14, 3)),
    ('Si IX', (14, 8)),
])
def test_species_tuple_to_species_string(species_string, species_tuple):
    assert species_tuple_to_string(species_tuple) == species_string
