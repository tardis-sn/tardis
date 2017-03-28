import os
import pytest
import  yaml
import h5py
import numpy as np

from tardis import __path__ as path
from tardis.io.util import yaml_load_config_file
from tardis.simulation import Simulation
from tardis.io.config_reader import Configuration
from tardis.plasma.base import *

@pytest.fixture(scope='module')
def partition_compare_data_fname():
    fname = 'partition_compare_data.h5'
    return os.path.join(path[0], 'plasma', 'tests', 'data', fname)

@pytest.fixture()
def partition_compare_data(partition_compare_data_fname):
    return h5py.File(partition_compare_data_fname, 'r')

@pytest.mark.skipif(not pytest.config.getvalue("partition-test"),
                    reason="partition tests are not included in this run")
@pytest.mark.partition
class TESTPartition():

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request, reference, data_path, pytestconfig):

        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
            pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(self.atom_data_filename), ("{0} atomic datafiles"
                                                         " does not seem to "
                                                         "exist".format(
            self.atom_data_filename))
        self.config_yaml = yaml_load_config_file('address')  #address
        self.config_yaml['atom_data'] = self.atom_data_filename
        conf = Configuration.from_config_dict(self.config_yaml)
        self.result = Simulation.from_config(conf,
                                             atom_data=self.atom_data).plasma    #error may come
        self.result.run()
        self.name = data_path['setup_name']

        if request.config.getoption("--generate-reference"):

            ref_data_path = os.path.join(
                data_path['gen_ref_path'], "{0}.h5".format(self.name)
            )
            if os.path.exists(ref_data_path):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(ref_data_path))
            self.result.to_hdf(path_or_buf=ref_data_path,
                               suffix_count=False)
            pytest.skip("Reference data saved at {0}".format(
                data_path['gen_ref_path']
            ))

    def test_level_boltzmann_factor_lte(level_boltzmann_factor_lte, partition_compare_data, levels):
        expected_lbfl = partition_compare_data['lbfl']
        assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[0], expected_lbfl[0])
        assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[1].ix[0], expected_lbfl[1])
        assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[10], expected_lbfl[2])

    def test_level_boltzmann_factor_dilute_lte(level_boltzmann_factor_dilute_lte, partition_compare_data,
        levels):
        expected_lbfdl = partition_compare_data['lbfdl']
        assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[0], expected_lbfdl[0])
        assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[1].ix[0], expected_lbfdl[1])
        assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[10], expected_lbfdl[2])

    def test_lte_partition_function(partition_function, partition_compare_data, levels):
        expected_pf = partition_compare_data['pf']
        assert np.allclose(partition_function.ix[2].ix[0], expected_pf[0])
        assert np.allclose(partition_function.ix[2].ix[1], expected_pf[1])
        assert np.allclose(partition_function.ix[2].ix[2], expected_pf[2])



