# tests for the config reader module

from tardis import config_reader, atomic
from astropy import units as u
import pytest



def test_quantity_parser_normal():
    q1 = config_reader.parse_quantity('5 km/s')
    assert q1.value == 5.
    assert q1.unit == u.Unit('km/s')

def test_quantity_parser_malformed_quantity1():
    with pytest.raises(config_reader.TardisConfigurationError):
        q1 = config_reader.parse_quantity('abcd')

def test_quantity_parser_malformed_quantity2():
    with pytest.raises(config_reader.TardisMalformedQuantityError):
        q1 = config_reader.parse_quantity('5 abcd')


def test_config_namespace_attribute_test():
    namespace = config_reader.TardisConfigurationNameSpace({'param1':1})
    assert namespace.param1 == 1

def test_config_namespace_attribute_test():
    namespace = config_reader.TardisConfigurationNameSpace({'param1':1})
    with pytest.raises(AttributeError):
        assert namespace.param2 == 1


def test_element_symbol_reformatter():
    def _test_element_symbol_reformatter(unformatted_element_string, formatted_element_string):
        assert config_reader.reformat_element_symbol(unformatted_element_string) == formatted_element_string

    data = [('si', 'Si'),
            ('sI', 'Si'),
            ('Si', 'Si'),
            ('c', 'C'),
            ('C', 'C'),
            ]

    for unformatted_element_string, formatted_element_string in data:
        yield _test_element_symbol_reformatter, unformatted_element_string, formatted_element_string


def test_element_symbol2atomic_number():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    def _test_element_symbol2atomic_number(element_string, atomic_number):
        assert config_reader.element_symbol2atomic_number(element_string, atom_data) == atomic_number

    data = [('sI', 14),
            ('ca', 20),
            ('Fe', 26)]

    for element_symbol, atomic_number in data:
        yield _test_element_symbol2atomic_number, element_symbol, atomic_number



def test_species_string_to_species():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    def _test_species_string_to_species_tuple(species_string, species_tuple):
        assert config_reader.species_string_to_tuple(species_string, atom_data) == species_tuple

    data = [('si ii', (14, 1) ),
            ('si 2', (14, 1)),
            ('si ix', (14, 8)),
            ]

    for species_string, species_tuple in data:
        yield _test_species_string_to_species_tuple, species_string, species_tuple
