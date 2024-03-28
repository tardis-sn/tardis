"""
Basic TARDIS Benchmark.
"""
from io import StringIO

import numpy as np
from astropy import units as u
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.util.base import (
    MalformedSpeciesError,
    MalformedElementSymbolError,
    MalformedQuantityError,
    int_to_roman,
    roman_to_int,
    calculate_luminosity,
    intensity_black_body,
    species_tuple_to_string,
    species_string_to_tuple,
    parse_quantity,
    element_symbol2atomic_number,
    atomic_number2element_symbol,
    reformat_element_symbol,
    quantity_linspace,
    convert_abundances_format,
)


# @skip_benchmark
class BenchmarkUtil(BenchmarkBase):
    """
    Class to benchmark the util function.
    """

    def __init__(self):
        pass

    @property
    def artis_abundances_fname(self):
        return f"{self.example_model_file_dir}/artis_abundances.dat"

    @staticmethod
    def time_malformed_species_error():
        malformed_species_error = MalformedSpeciesError("He")
        assert malformed_species_error.malformed_element_symbol == "He"
        assert (
                str(malformed_species_error)
                == 'Expecting a species notation (e.g. "Si 2", "Si II", "Fe IV") - supplied He'
        )

    @staticmethod
    def time_malformed_elements_symbol_error():
        malformed_elements_symbol_error = MalformedElementSymbolError("Hx")
        assert malformed_elements_symbol_error.malformed_element_symbol == "Hx"
        assert (
                str(malformed_elements_symbol_error)
                == "Expecting an atomic symbol (e.g. Fe) - supplied Hx"
        )

    @staticmethod
    def time_malformed_quantity_error():
        malformed_quantity_error = MalformedQuantityError("abcd")
        assert malformed_quantity_error.malformed_quantity_string == "abcd"
        assert (
                str(malformed_quantity_error)
                == 'Expecting a quantity string(e.g. "5 km/s") for keyword - supplied abcd'
        )

    @staticmethod
    @parameterize({"Parameters": [
        {
            "value": 1,
            "expected": "I"
        },
        {
            "value": 5,
            "expected": "V"
        },
        {
            "value": 19,
            "expected": "XIX"
        },
        {
            "value": 556,
            "expected": "DLVI"
        },
        {
            "value": 1400,
            "expected": "MCD"
        },
        {
            "value": 1999,
            "expected": "MCMXCIX"
        },
        {
            "value": 3000,
            "expected": "MMM"
        },
    ]})
    def time_int_to_roman(parameters):
        value = parameters["value"]
        expected = parameters["expected"]
        assert int_to_roman(value) == expected

    @staticmethod
    @parameterize({"Parameters": [
        {
            "value": "I",
            "expected": 1
        },
        {
            "value": "V",
            "expected": 5
        },
        {
            "value": "XIX",
            "expected": 19
        },
        {
            "value": "DLVI",
            "expected": 556
        },
        {
            "value": "MCD",
            "expected": 1400
        },
        {
            "value": "MCMXCIX",
            "expected": 1999
        },
        {
            "value": "MMM",
            "expected": 3000
        },
    ]})
    def time_roman_to_int(parameters):
        value = parameters["value"]
        expected = parameters["expected"]

        assert roman_to_int(value) == expected

    @staticmethod
    @skip_benchmark
    @parameterize({"Parameters": [
        {
            "string_io": StringIO("4000 1e-21\n4500 3e-21\n5000 5e-21"),
            "distance": "100 km",
            "result": (0.0037699111843077517, 4000.0, 5000.0)
        },
        {
            "string_io": StringIO("7600 2.4e-19\n7800 1.6e-19\n8100 9.1e-20"),
            "distance": "500 km",
            "result": (2.439446695512474, 7600.0, 8100.0)
        },
    ]})
    def time_calculate_luminosity(parameters):
        string_io = parameters["string_io"]
        distance = parameters["distance"]
        result = parameters["result"]

        assert calculate_luminosity(string_io, distance) == result

    @staticmethod
    @parameterize({"Parameters": [
        {
            "nu": 10 ** 6,
            "t": 1000,
            "i": 3.072357852080765e-22
        },
        {
            "nu": 10 ** 6,
            "t": 300,
            "i": 9.21707305730458e-23
        },
        {
            "nu": 10 ** 8,
            "t": 1000,
            "i": 6.1562660718558254e-24
        },
        {
            "nu": 10 ** 8,
            "t": 300,
            "i": 1.846869480674048e-24
        },
    ]})
    def time_intensity_black_body(parameters):
        nu = parameters["nu"]
        t = parameters["t"]
        i = parameters["i"]
        assert np.isclose(intensity_black_body(nu, t), i)

    @staticmethod
    @parameterize({"Parameters": [
        {
            "species_tuple": (14, 1),
            "roman_numerals": True,
            "species_string": "Si II"
        },
        {
            "species_tuple": (14, 1),
            "roman_numerals": False,
            "species_string": "Si 1"
        },
        {
            "species_tuple": (14, 3),
            "roman_numerals": True,
            "species_string": "Si IV"
        },
        {
            "species_tuple": (14, 3),
            "roman_numerals": False,
            "species_string": "Si 3"
        },
        {
            "species_tuple": (14, 8),
            "roman_numerals": True,
            "species_string": "Si IX"
        },
        {
            "species_tuple": (14, 8),
            "roman_numerals": False,
            "species_string": "Si 8"
        },
    ]})
    def time_species_tuple_to_string(parameters):
        species_tuple = parameters["species_tuple"]
        roman_numerals = parameters["roman_numerals"]
        species_string = parameters["species_string"]
        assert (
                species_tuple_to_string(species_tuple, roman_numerals=roman_numerals)
                == species_string
        )

    @staticmethod
    @parameterize({"Parameters": [
        {
            "species_string": "si ii",
            "species_tuple": (14, 1)
        },
        {
            "species_string": "si 2",
            "species_tuple": (14, 1)
        },
        {
            "species_string": "si ix",
            "species_tuple": (14, 8)
        },
    ]})
    def time_species_string_to_tuple(parameters):
        species_string = parameters["species_string"]
        species_tuple = parameters["species_tuple"]
        assert species_string_to_tuple(species_string) == species_tuple

        try:
            species_string_to_tuple("II")
            assert False
        except MalformedSpeciesError:
            assert True
        except:
            assert False

        try:
            species_string_to_tuple("He Si")
            assert False
        except MalformedSpeciesError:
            assert True
        except:
            assert False

        try:
            species_string_to_tuple("He IX")
            assert False
        except ValueError:
            assert True
        except:
            assert False

    @staticmethod
    def time_parse_quantity():
        q1 = parse_quantity("5 km/s")
        assert q1.value == 5.0
        assert q1.unit == u.Unit("km/s")

        try:
            parse_quantity(5)
            assert False
        except MalformedQuantityError:
            assert True
        except:
            assert False

        try:
            parse_quantity("abcd")
            assert False
        except MalformedQuantityError:
            assert True
        except:
            assert False

        try:
            parse_quantity("a abcd")
            assert False
        except MalformedQuantityError:
            assert True
        except:
            assert False

        try:
            parse_quantity("5 abcd")
            assert False
        except MalformedQuantityError:
            assert True
        except:
            assert False

    @staticmethod
    @parameterize({"Parameters": [
        {
            "element_symbol": "sI",
            "atomic_number": 14
        },
        {
            "element_symbol": "ca",
            "atomic_number": 20
        },
        {
            "element_symbol": "Fe",
            "atomic_number": 26
        },
    ]})
    def time_element_symbol2atomic_number(parameters):
        element_symbol = parameters["element_symbol"]
        atomic_number = parameters["atomic_number"]
        assert element_symbol2atomic_number(element_symbol) == atomic_number

        try:
            element_symbol2atomic_number("Hx")
            assert False
        except MalformedElementSymbolError:
            assert True
        except:
            assert False

    @staticmethod
    def time_atomic_number2element_symbol():
        assert atomic_number2element_symbol(14) == "Si"

    @staticmethod
    @parameterize({"Parameters": [
        {
            "unformatted_element_string": "si",
            "formatted_element_string": "Si"
        },
        {
            "unformatted_element_string": "sI",
            "formatted_element_string": "Si"
        },
        {
            "unformatted_element_string": "Si",
            "formatted_element_string": "Si"
        },
        {
            "unformatted_element_string": "c",
            "formatted_element_string": "C"
        },
        {
            "unformatted_element_string": "C",
            "formatted_element_string": "C"
        },
    ]})
    def time_reformat_element_symbol(parameters):
        unformatted_element_string = parameters["unformatted_element_string"]
        formatted_element_string = parameters["formatted_element_string"]
        assert (
                reformat_element_symbol(unformatted_element_string)
                == formatted_element_string
        )

    @staticmethod
    @parameterize({"Parameters": [
        {
            "start": u.Quantity(1, "km/s"),
            "stop": u.Quantity(5, "km/s"),
            "num": 5,
            "expected": u.Quantity(np.array([1.0, 2.0, 3.0, 4.0, 5.0]), "km/s")
        },
        {
            "start": u.Quantity(0.5, "eV"),
            "stop": u.Quantity(0.6, "eV"),
            "num": 3,
            "expected": u.Quantity(np.array([0.5, 0.55, 0.6]), "eV")
        },
    ]})
    def time_quantity_linspace(parameters):
        start = parameters["start"]
        stop = parameters["stop"]
        num = parameters["num"]
        expected = parameters["expected"]
        obtained = quantity_linspace(start, stop, num)
        assert obtained.unit == expected.unit
        assert obtained.value.all() == expected.value.all()

        try:
            quantity_linspace(u.Quantity(0.5, "eV"), "0.6 eV", 3)
            assert False
        except ValueError:
            assert True
        except:
            assert False

    def time_convert_abundances_format(self):
        artis_abundances_fname = self.artis_abundances_fname
        abundances = convert_abundances_format(artis_abundances_fname)
        assert np.isclose(abundances.loc[3, "O"], 1.240199e-08, atol=1.0e-12)
        assert np.isclose(abundances.loc[1, "Co"], 2.306023e-05, atol=1.0e-12)
        assert np.isclose(abundances.loc[69, "Ni"], 1.029928e-17, atol=1.0e-12)
        assert np.isclose(abundances.loc[2, "C"], 4.425876e-09, atol=1.0e-12)
