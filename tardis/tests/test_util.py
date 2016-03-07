import pytest
import numpy as np
from astropy import units as u
from io import StringIO

from tardis.util import MalformedSpeciesError, MalformedElementSymbolError, MalformedQuantityError
from tardis.util import (int_to_roman, roman_to_int, create_synpp_yaml,
                         calculate_luminosity, intensity_black_body, savitzky_golay,
                         species_tuple_to_string, species_string_to_tuple, parse_quantity,
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


@pytest.mark.parametrize(['test_input', 'expected_result'], [
    (1, 'I'),
    (5, 'V'),
    (19, 'XIX'),
    (556, 'DLVI'),
    (1400, 'MCD'),
    (1999, 'MCMXCIX'),
    (3000, 'MMM')
])
def test_int_to_roman(test_input, expected_result):
    assert int_to_roman(test_input) == expected_result

    with pytest.raises(TypeError):
        int_to_roman(1.5)

    with pytest.raises(ValueError):
        int_to_roman(0)

    with pytest.raises(ValueError):
        int_to_roman(4000)


@pytest.mark.parametrize(['test_input', 'expected_result'], [
    ('I', 1),
    ('V', 5),
    ('XIX', 19),
    ('DLVI', 556),
    ('MCD', 1400),
    ('MCMXCIX', 1999),
    ('MMM', 3000)
])
def test_roman_to_int(test_input, expected_result):
    assert roman_to_int(test_input) == expected_result

    with pytest.raises(TypeError):
        roman_to_int(1)

    with pytest.raises(ValueError):
        roman_to_int('IIV')


@pytest.mark.parametrize(['string_io', 'distance', 'result'], [
    (StringIO(u'4000 1e-21\n4500 3e-21\n5000 5e-21'), '100 km', (0.0037699111843077517, 4000.0, 5000.0)),
    (StringIO(u'7600 2.4e-19\n7800 1.6e-19\n8100 9.1e-20'), '500 km', (2.439446695512474, 7600.0, 8100.0))
])
def test_calculate_luminosity(string_io, distance, result):
    assert calculate_luminosity(string_io, distance) == result


@pytest.mark.parametrize(['nu', 't', 'i'], [
    (10**6, 1000, 3.072357852080765e-22),
    (10**6, 300, 9.21707305730458e-23),
    (10**8, 1000, 6.1562660718558254e-24),
    (10**8, 300, 1.846869480674048e-24),
])
def test_intensity_black_body(nu, t, i):
    assert float(intensity_black_body(nu, t)) == i


def test_savitzky_golay():
    # simple testcase:
    # time axis sampled into 5 values between 0 and 0.5
    t = np.linspace(0, 0.5, 5)
    # y axis represent signal values - time complexed
    # provided undisturbed signal is ramp input, some disturbances are added randomly
    y = np.exp(t) + np.array([0.00016543, -0.00011681, -0.00060518, -0.00020232, 0.0006262])

    # applying the filter
    ysg = savitzky_golay(y, window_size=21, order=5)

    # expected result on application of filter
    result = np.array([0.62999136, 0.6976977, 0.93429473, 1.23388831, 1.51577056,
                       1.72441978, 1.82950054, 1.82586363, 1.73354608])
    assert ysg.all() == result.all()

    with pytest.raises(ValueError):
        # window_size and order have to be of type int
        savitzky_golay(y, window_size='a', order='b')

    with pytest.raises(TypeError):
        # window_size size must be a positive odd number
        savitzky_golay(y, window_size=10, order=4)

    with pytest.raises(TypeError):
        # window_size is too small for the polynomials order
        savitzky_golay(y, window_size=1, order=4)


@pytest.mark.parametrize(['species_tuple', 'roman_numerals', 'species_string'], [
    ((14, 1), True, 'Si II'),
    ((14, 1), False, 'Si 1'),
    ((14, 3), True, 'Si IV'),
    ((14, 3), False, 'Si 3'),
    ((14, 8), True, 'Si IX'),
    ((14, 8), False, 'Si 8'),
])
def test_species_tuple_to_string(species_tuple, roman_numerals, species_string):
    assert species_tuple_to_string(species_tuple, roman_numerals=roman_numerals) == species_string


@pytest.mark.parametrize(['species_string', 'species_tuple'], [
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


@pytest.mark.parametrize(['element_symbol', 'atomic_number'], [
    ('sI', 14),
    ('ca', 20),
    ('Fe', 26)
])
def test_element_symbol2atomic_number(element_symbol, atomic_number):
    assert element_symbol2atomic_number(element_symbol) == atomic_number

    with pytest.raises(MalformedElementSymbolError):
        element_symbol2atomic_number('Hx')


def test_atomic_number2element_symbol():
    assert atomic_number2element_symbol(14) == 'Si'


@pytest.mark.parametrize(['unformatted_element_string', 'formatted_element_string'], [
    ('si', 'Si'),
    ('sI', 'Si'),
    ('Si', 'Si'),
    ('c', 'C'),
    ('C', 'C'),
])
def test_reformat_element_symbol(unformatted_element_string, formatted_element_string):
    assert reformat_element_symbol(unformatted_element_string) == formatted_element_string


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
