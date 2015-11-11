import pytest
import numpy as np
import yaml
import tardis
import numpy.testing as nptesting
from astropy import units as u
import os

from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel

from tardis.io.config_reader import Configuration



def data_path(fname):
    return os.path.join(tardis.__path__[0], 'tests', 'data', fname)

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
class TestSimpleRun():
    """
    Very simple run
    """

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
            'tardis/io/tests/data/tardis_configv1_verysimple.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        tardis_config = Configuration.from_config_dict(self.config_yaml)
        self.model = Radial1DModel(tardis_config)
        self.simulation = Simulation(tardis_config)

        self.simulation.legacy_run_simulation(self.model)

    def test_spectrum(self):
        luminosity_density = np.load(
            data_path('simple_test_luminosity_density_lambda.npy'))

        luminosity_density = luminosity_density * u.Unit(
            'erg / (Angstrom s)')

        np.testing.assert_allclose(
            self.model.spectrum.luminosity_density_lambda,luminosity_density)
