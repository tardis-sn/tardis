import pytest
import numpy as np
import yaml
import tardis
from tardis import io, model, simulation
from tardis.io.config_reader import Configuration
from numpy.testing import assert_array_almost_equal
from astropy import units as u
import os


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
        assert os.path.exists(self.atom_data_filename), ("{0} atomic datafiles "
                                                         "does not seem to "
                                                         "exist".format(
            self.atom_data_filename))
        self.config_yaml = yaml.load(open('tardis/io/tests/data/tardis_configv1_verysimple.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        self.config = Configuration.from_config_dict(self.config_yaml)
        self.model = model.Radial1DModel(self.config)
        simulation.run_radial1d(self.model)


    def test_spectrum(self):
        wavelength, luminosity_density = \
            np.loadtxt(data_path('simple_test_spectrum.dat'), unpack=True)

        luminosity_density = luminosity_density * u.Unit('erg / (Angstrom s)')

        assert_array_almost_equal(self.model.spectrum.luminosity_density_lambda,
                                  luminosity_density)
