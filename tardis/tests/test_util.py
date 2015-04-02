#tests for util module

import pytest
from astropy import units as u
from tardis import atomic
from tardis.util import species_string_to_tuple, parse_quantity, element_symbol2atomic_number, atomic_number2element_symbol, reformat_element_symbol, MalformedQuantityError, roman_to_int, int_to_roman

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

def test_element_symbol2atomic_number():
    def _test_element_symbol2atomic_number(element_string, atomic_number):
        assert element_symbol2atomic_number(element_string) == atomic_number

    data = [('sI', 14),
            ('ca', 20),
            ('Fe', 26)]

    for element_symbol, atomic_number in data:
        yield _test_element_symbol2atomic_number, element_symbol, atomic_number

def test_element_symbol_reformatter():
    def _test_element_symbol_reformatter(unformatted_element_string, formatted_element_string):
        assert reformat_element_symbol(unformatted_element_string) == formatted_element_string

    data = [('si', 'Si'),
            ('sI', 'Si'),
            ('Si', 'Si'),
            ('c', 'C'),
            ('C', 'C'),
            ]

    for unformatted_element_string, formatted_element_string in data:
        yield _test_element_symbol_reformatter, unformatted_element_string, formatted_element_string

def test_species_string_to_species():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    def _test_species_string_to_species_tuple(species_string, species_tuple):
        assert species_string_to_tuple(species_string) == species_tuple

    data = [('si ii', (14, 1) ),
            ('si 2', (14, 1)),
            ('si ix', (14, 8)),
            ]

    for species_string, species_tuple in data:
        yield _test_species_string_to_species_tuple, species_string, species_tuple

testData = [('CMXC', 990),
	    ('CCCLVI', 356),
	    ('MCCCXV', 1315), 
	    ('DCLXXIV', 674),
	   ]

@pytest.mark.parametrize("roman_number, int_number", testData)
def test_roman_to_int(roman_number, int_number):
    assert roman_to_int(roman_number) == int_number

testData = [(823, 'DCCCXXIII'),
	    (2421, 'MMCDXXI'),
	    (372, 'CCCLXXII'),
	    (1443, 'MCDXLIII'),
	   ]

@pytest.mark.parametrize("int_number, roman_number", testData)
def test_int_to_roman(int_number, roman_number):
    assert int_to_roman(int_number) == roman_number

"""This test is written to check the boundary conditions for the functions
   roman_to_int and int_to_roman"""
def test_boundary():
    with pytest.raises(ValueError):
	int_number = roman_to_int("MMMM")

    with pytest.raises(ValueError):
	roman_number = int_to_roman("4000")
