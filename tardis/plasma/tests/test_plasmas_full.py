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
            'tardis/plasma/tests/plasma_test_config_lte.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename
        self.lte_model = run_tardis(self.config_yaml)

    def test_lte_plasma(self):
        old_plasma_t_rads = \
            np.loadtxt(data_path('plasma_comparison_lte_trads.dat'),
                unpack=True)

        old_plasma_t_rads = old_plasma_t_rads * u.Unit('K')

        np.testing.assert_allclose(
            self.lte_model.t_rads, old_plasma_t_rads, atol=50 * u.Unit('K'))
