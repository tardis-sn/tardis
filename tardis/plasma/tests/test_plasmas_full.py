import pytest
import numpy as np
import yaml
import tardis
import numpy.testing as nptesting
from astropy import units as u
import os

from tardis.base import run_tardis

def data_path(fname):
    return os.path.join(tardis.__path__[0], 'plasma', 'tests', 'data', fname)

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
        self.config_yaml = yaml.load(open(
            'tardis/plasma/tests/data/plasma_test_config_lte.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename
        self.lte_model = run_tardis(self.config_yaml)
        self.config_yaml = yaml.load(open(
            'tardis/plasma/tests/data/plasma_test_config_nlte.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename
        self.nlte_model = run_tardis(self.config_yaml)

    def test_lte_plasma(self):
        old_plasma_t_rads = \
            np.loadtxt(data_path('plasma_comparison_lte_trads.dat'),
                unpack=True)
        new_plasma_t_rads = self.lte_model.t_rads / u.Unit('K')
        old_plasma_levels = \
            np.loadtxt(data_path('plasma_comparison_lte_levels.dat'),
                unpack=True)
        new_plasma_levels = \
            self.lte_model.plasma_array.get_value(
            'level_number_density').ix[8].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=100)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.05)

    def test_nlte_plasma(self):
        old_plasma_t_rads = \
            np.loadtxt(data_path('plasma_comparison_nlte_trads.dat'),
                unpack=True)
        new_plasma_t_rads = self.nlte_model.t_rads / u.Unit('K')
        old_plasma_levels = \
            np.loadtxt(data_path('plasma_comparison_nlte_levels.dat'),
                unpack=True)
        new_plasma_levels = \
            self.nlte_model.plasma_array.get_value(
            'level_number_density').ix[2].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=150)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1)