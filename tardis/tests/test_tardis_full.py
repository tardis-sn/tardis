import pytest
import numpy as np
import yaml
import tardis
from tardis import io, model, simulation
from tardis.io.config_reader import Configuration
import numpy.testing as nptesting
from astropy import units as u
import os


from tardis.base import run_tardis

def data_path(fname):
    return os.path.join(tardis.__path__[0], 'tests', 'data', fname)

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
@pytest.mark.xfail
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

        self.model = run_tardis(self.config_yaml)


    def test_spectrum(self):
        wavelength, luminosity_density = \
            np.loadtxt(data_path('simple_test_spectrum.dat'), unpack=True)

        luminosity_density = luminosity_density * u.Unit('erg / (Angstrom s)')

        np.testing.assert_allclose(self.model.spectrum.luminosity_density_lambda,
                                  luminosity_density)
