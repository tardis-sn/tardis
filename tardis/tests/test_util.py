#tests for util module

import pytest
from astropy import units as u
from tardis import atomic
from tardis.util import species_string_to_tuple, parse_quantity, element_symbol2atomic_number, atomic_number2element_symbol, reformat_element_symbol, MalformedQuantityError, species_tuple_to_string

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

testData = [('sI', 14),
        ('ca', 20),
        ('Fe', 26)]

@pytest.mark.parametrize("element_string, atomic_number", testData)
def test_element_symbol2atomic_number(element_string, atomic_number):
    assert element_symbol2atomic_number(element_string) == atomic_number

testData = [('si', 'Si'),
            ('sI', 'Si'),
            ('Si', 'Si'),
            ('c', 'C'),
            ('C', 'C'),
           ]

@pytest.mark.parametrize("unformatted_element_string, formatted_element_string", testData)
def test_element_symbol_reformatter(unformatted_element_string, formatted_element_string):
    assert reformat_element_symbol(unformatted_element_string) == formatted_element_string

testData = [('si ii', (14, 1) ),
            ('si 2', (14, 1)),
            ('si ix', (14, 8)),
           ]

@pytest.mark.parametrize("species_string, species_tuple", testData)
def test_species_string_to_species_tuple(species_string, species_tuple):
    assert species_string_to_tuple(species_string) == species_tuple

testData = [((14, 9), True, 'Si X'),
            ((14, 100), False, 'Si 100'),
	    ((14, 100), True, 'Si CI'),
           ]

@pytest.mark.parametrize("species_tuple, roman_numerals, species_string", testData)
def test_species_tuple_to_string(species_tuple, roman_numerals, species_string):
    assert species_tuple_to_string(species_tuple, roman_numerals) == species_string
