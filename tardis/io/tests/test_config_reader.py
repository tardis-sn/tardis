# tests for the config reader module

from tardis import atomic
from tardis.io import config_reader
from astropy import units as u
import os
import pytest
import yaml

from numpy.testing import assert_almost_equal, assert_array_almost_equal
def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, 'data', filename)

"""
These two are the same function name, but do havev a slight difference
"""
def test_config_namespace_attribute_test():
    namespace = config_reader.TARDISConfigurationNameSpace({'param1':1})
    assert namespace.param1 == 1

def test_config_namespace_attribute_test():
    namespace = config_reader.TARDISConfigurationNameSpace({'param1':1})
    with pytest.raises(AttributeError):
        assert namespace.param2 == 1

"""
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
"""

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
        assert_almost_equal(self.config.abundances.ix[14, 31], 2.156751e-01)

     
class TestParseConfigV1ArtisDensityAbundancesAllAscii:

    def setup(self):
        #general parsing of the paper config
        filename = 'tardis_configv1_ascii_density_abund.yml'

        self.yaml_data = yaml.load(open(data_path(filename)))
        self.yaml_data['model']['structure']['filename'] = 'tardis/io/tests/data/density.dat'
        self.yaml_data['model']['abundances']['filename'] = 'tardis/io/tests/data/abund.dat'
    
        self.config = config_reader.TARDISConfiguration.from_config_dict(self.yaml_data,
                                                                  test_parser=True)


    def test_velocities(self):
        assert self.config.structure.v_inner.unit == u.Unit('cm/s')
        assert_almost_equal(self.config.structure.v_inner[0], 11000 * u.km/u.s)

    def test_abundances(self):
        assert_almost_equal(self.config.abundances.ix[14, 0], 0.1)
        assert_almost_equal(self.config.abundances.ix[14, 1], 0.2)
        assert_almost_equal(self.config.abundances.ix[14, 2], 0.2)
        assert_almost_equal(self.config.abundances.ix[14, 3], 0.2)
        assert_almost_equal(self.config.abundances.ix[14, 4], 0.2)
        assert_almost_equal(self.config.abundances.ix[14, 5], 0.2)
        assert_almost_equal(self.config.abundances.ix[14, 6], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 0], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 1], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 2], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 3], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 4], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 5], 0.0)
        assert_almost_equal(self.config.abundances.ix[6, 6], 0.5)

    def test_densities(self):
        assert_almost_equal(self.config.structure.mean_densities[0], 9.7656229e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[1], 4.8170911e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[2], 2.5600000e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[3], 1.4450533e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[4], 8.5733893e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[5], 5.3037103e-11 * u.Unit('g/cm3') / 13.0**3 )
        assert_almost_equal(self.config.structure.mean_densities[6], 3.3999447e-11 * u.Unit('g/cm3') / 13.0**3 )

#write tests for inner and outer boundary indices
