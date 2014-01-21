# tests for the config reader module

from tardis import atomic
from tardis.io import config_reader
from astropy import units as u
import os
import pytest
import yaml

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from tardis.util import species_string_to_tuple, parse_quantity, element_symbol2atomic_number, reformat_element_symbol, MalformedQuantityError


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, 'data', filename)



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


def test_config_namespace_attribute_test():
    namespace = config_reader.TARDISConfigurationNameSpace({'param1':1})
    assert namespace.param1 == 1

def test_config_namespace_attribute_test():
    namespace = config_reader.TARDISConfigurationNameSpace({'param1':1})
    with pytest.raises(AttributeError):
        assert namespace.param2 == 1


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


def test_element_symbol2atomic_number():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    def _test_element_symbol2atomic_number(element_string, atomic_number):
        assert element_symbol2atomic_number(element_string) == atomic_number

    data = [('sI', 14),
            ('ca', 20),
            ('Fe', 26)]

    for element_symbol, atomic_number in data:
        yield _test_element_symbol2atomic_number, element_symbol, atomic_number

def test_quantity_linspace():
    quantity_linspace_dict = dict(start='1.1e4 km/s', stop='2e4 cm/h', num=1000)
    quantity_linspace = config_reader.parse_quantity_linspace(quantity_linspace_dict)
    assert_almost_equal(quantity_linspace[0].value, 1.1e4)
    assert_almost_equal(quantity_linspace[-1].to('cm/h').value, 2e4)
    assert len(quantity_linspace) == 1001


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

class TestParsePaper1Config:

    def setup(self):
        #general parsing of the paper config
        self.config = config_reader.TARDISConfiguration.from_yaml(data_path('paper1_tardis_configv1.yml'),
                                                                  test_parser=True)
        self.yaml_data = yaml.load(open(data_path('paper1_tardis_configv1.yml')))



    def test_abundances(self):
        oxygen_abundance = self.yaml_data['model']['abundances']['O']
        assert_array_almost_equal(oxygen_abundance, self.config.abundances.ix[8].values)

        assert True

    def test_velocities(self):
        assert_almost_equal(parse_quantity(self.yaml_data['model']['structure']['velocity']['start']),
                            self.config.structure.v_inner[0])
        assert_almost_equal(parse_quantity(self.yaml_data['model']['structure']['velocity']['stop']),
                    self.config.structure.v_outer[-1])
        assert len(self.config.structure.v_outer) == self.yaml_data['model']['structure']['velocity']['num']

    def test_densities(self):
        pass

class TestParseConfigV1ASCIIDensity:

    def setup(self):
        #general parsing of the paper config
        filename = 'tardis_configv1_ascii_density.yml'
        self.config = config_reader.TARDISConfiguration.from_yaml(data_path(filename),
                                                                  test_parser=True)
        self.yaml_data = yaml.load(open(data_path(filename)))


    def test_velocities(self):
        assert self.config.structure.v_inner.unit == u.Unit('cm/s')
        assert_almost_equal(self.config.structure.v_inner[0].value, 1e4 * 1e5)

    def test_abundances(self):
        oxygen_abundance = self.yaml_data['model']['abundances']['O']
        assert_array_almost_equal(oxygen_abundance, self.config.abundances.ix[8].values)


class TestParseConfigV1ArtisDensity:

    def setup(self):
        #general parsing of the paper config
        filename = 'tardis_configv1_artis_density.yml'
        self.config = config_reader.TARDISConfiguration.from_yaml(data_path(filename),
                                                                  test_parser=True)
        self.yaml_data = yaml.load(open(data_path(filename)))


    def test_velocities(self):
        assert self.config.structure.v_inner.unit == u.Unit('cm/s')
        assert_almost_equal(self.config.structure.v_inner[0].value, 1.259375e+03 * 1e5)

    def test_abundances(self):
        oxygen_abundance = self.yaml_data['model']['abundances']['O']
        assert_array_almost_equal(oxygen_abundance, self.config.abundances.ix[8].values)




class TestParseConfigV1ArtisDensityAbundances:

    def setup(self):
        #general parsing of the paper config
        filename = 'tardis_configv1_artis_density.yml'

        self.yaml_data = yaml.load(open(data_path(filename)))
        self.yaml_data['model']['abundances'] = {'type': 'file',
                                                 'filename': 'tardis/io/tests/data/artis_abundances.dat',
                                                 'filetype': 'artis'}

        self.config = config_reader.TARDISConfiguration.from_config_dict(self.yaml_data,
                                                                  test_parser=True)


    def test_velocities(self):
        assert self.config.structure.v_inner.unit == u.Unit('cm/s')
        assert_almost_equal(self.config.structure.v_inner[0].value, 1.259375e+03 * 1e5)

    def test_abundances(self):
        assert_almost_equal(self.config.abundances.ix[14, 54], 0.21864420000000001)


class TestParseConfigV1ArtisDensityAbundancesVSlice:

    def setup(self):
        #general parsing of the paper config
        filename = 'tardis_configv1_artis_density_v_slice.yml'

        self.yaml_data = yaml.load(open(data_path(filename)))
        self.yaml_data['model']['abundances'] = {'type': 'file',
                                                 'filename': 'tardis/io/tests/data/artis_abundances.dat',
                                                 'filetype': 'artis'}

        self.config = config_reader.TARDISConfiguration.from_config_dict(self.yaml_data,
                                                                  test_parser=True)


    def test_velocities(self):
        assert self.config.structure.v_inner.unit == u.Unit('cm/s')
        assert_almost_equal(self.config.structure.v_inner[0], 9000 * u.km/u.s)

    def test_abundances(self):
        assert_almost_equal(self.config.abundances.ix[14, 54], 0.21864420000000001)



#write tests for inner and outer boundary indices