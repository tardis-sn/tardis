import pytest
import numpy as np
import tardis
import numpy.testing as nptesting
from astropy import units as u
import os
import h5py

from tardis.base import run_tardis
from tardis.io.util import yaml_load_config_file


def data_path(fname):
    return os.path.join(tardis.__path__[0], 'plasma', 'tests', 'data', fname)

@pytest.fixture()
def plasma_compare_data_fname():
    return data_path('plasma_test_data.h5')

@pytest.fixture()
def plasma_compare_data(plasma_compare_data_fname):
    return h5py.File(plasma_compare_data_fname, 'r')

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
class TestPlasmas():
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
            pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(self.atom_data_filename), ("{0} atomic datafiles"
                                                         " does not seem to "
                                                         "exist".format(
            self.atom_data_filename))
        self.config_yaml = yaml_load_config_file(
            'tardis/plasma/tests/data/plasma_test_config_lte.yml')
        self.config_yaml['atom_data'] = self.atom_data_filename
        self.lte_model = run_tardis(self.config_yaml)
        self.config_yaml = yaml_load_config_file(
            'tardis/plasma/tests/data/plasma_test_config_nlte.yml')
        self.config_yaml['atom_data'] = self.atom_data_filename
        self.nlte_model = run_tardis(self.config_yaml)

    def test_lte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data['test_lte1/t_rad']
        old_plasma_levels = plasma_compare_data['test_lte1/levels']

        new_plasma_t_rads = self.lte_model.t_rads / u.Unit('K')
        new_plasma_levels = \
            self.lte_model.plasma.get_value(
            'level_number_density').ix[8].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=100)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1)

    def test_nlte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data['test_nlte1/t_rad']
        old_plasma_levels = plasma_compare_data['test_nlte1/levels']
        new_plasma_t_rads = self.nlte_model.t_rads / u.Unit('K')
        new_plasma_levels = \
            self.nlte_model.plasma.get_value(
            'level_number_density').ix[2].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=150)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1)
