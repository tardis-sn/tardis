import pytest
import numpy as np
import tardis
import numpy.testing as nptesting
from astropy import units as u
import os

from tardis.io.util import yaml_load_config_file
from tardis.simulation.base import Simulation
from tardis.io.config_reader import Configuration



def data_path(fname):
    return os.path.join(tardis.__path__[0], 'tests', 'data', fname)

@pytest.mark.skipif(not pytest.config.getvalue("tardis-refdata"),
                    reason='--tardis-refdata was not specified')
class TestSimpleRun():
    """
    Very simple run
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
            os.path.join(pytest.config.getvalue('tardis-refdata'), 'atom_data', 'kurucz_cd23_chianti_H_He.h5')))
        assert os.path.exists(self.atom_data_filename), ("{0} atomic datafiles"
                                                         " does not seem to "
                                                         "exist".format(
            self.atom_data_filename))
        self.config_yaml = yaml_load_config_file(
            'tardis/io/tests/data/tardis_configv1_verysimple.yml')
        self.config_yaml['atom_data'] = self.atom_data_filename

        tardis_config = Configuration.from_config_dict(self.config_yaml)
        self.simulation = Simulation.from_config(tardis_config)
        self.simulation.run()

    def test_j_blue_estimators(self):
        j_blue_estimator = np.load(
            data_path('simple_test_j_blue_estimator.npy'))

        np.testing.assert_allclose(self.simulation.runner.j_blue_estimator,
                                   j_blue_estimator)

    def test_spectrum(self):
        luminosity_density = np.load(
            data_path('simple_test_luminosity_density_lambda.npy'))

        luminosity_density = luminosity_density * u.Unit(
            'erg / (Angstrom s)')

        np.testing.assert_allclose(
            self.simulation.runner.spectrum.luminosity_density_lambda,
            luminosity_density)

    def test_virtual_spectrum(self):
        virtual_luminosity_density = np.load(
            data_path('simple_test_virtual_luminosity_density_lambda.npy'))

        virtual_luminosity_density = virtual_luminosity_density * u.Unit(
            'erg / (Angstrom s)')

        np.testing.assert_allclose(
            self.simulation.runner.spectrum_virtual.luminosity_density_lambda,
            virtual_luminosity_density)

    def test_plasma_properties(self):

        pass

    def test_runner_properties(self):
        """Tests whether a number of runner attributes exist and also verifies
        their types

        Currently, runner attributes needed to call the model routine to_hdf5
        are checked.

        """

        virt_type = np.ndarray


        props_required_by_modeltohdf5 = dict([
                ("virt_packet_last_interaction_type", virt_type),
                ("virt_packet_last_line_interaction_in_id", virt_type),
                ("virt_packet_last_line_interaction_out_id", virt_type),
                ("virt_packet_last_interaction_in_nu", virt_type),
                ("virt_packet_nus", virt_type),
                ("virt_packet_energies", virt_type),
                ])

        required_props = props_required_by_modeltohdf5.copy()

        for prop, prop_type in required_props.items():

            assert type(getattr(self.simulation.runner, prop)) == prop_type, ("wrong type of attribute '{}': expected {}, found {}".format(prop, prop_type, type(getattr(self.simulation.runner, prop))))
