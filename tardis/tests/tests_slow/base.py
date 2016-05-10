import os
import yaml
import numpy as np
import numpy.testing as nptesting
from astropy import units as u

from collector import data_path
from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


class SlowTest(object):
    """
    Slow integration test for various setups under $TARDIS_PATH/tests/tests_slow.
    """

    @classmethod
    def base_setup(self, model_dir, config_file, abundances=None, densities=None):
        """
        This method does an initial setup of making the configuration and performing
        a single run of integration test.

        Args:
            model_dir: (~str)
                Particular directory name in `tests_slow` directory containing files
                for a specific setup.
            config_file: (~str or ~dict)
                Name of yml config file inside model_dir, or a dict containing config.
            abundances: (~str)
                Name of the data file containing abundances profile.
            densities: (~str)
                Name of the data file containing densities profile.
        """
        self.model_dir = model_dir
        self.config_file = config_file
        self.abundances = abundances
        self.densities = densities

        # First we check whether the config file does even exist at the desired path.
        assert os.path.exists(data_path('w7_13d', 'tardis_w7.yml')), \
            ('%s config file does not exist' % data_path(self.model_dir, self.config_file))

        # The available config file doesn't have the file paths of atom data file,
        # densities and abundances profile files as desired. We form a dictionary
        # from the config file and override those parameters - putting file paths
        # of these three files at proper places.
        config_yaml = yaml.load(open(data_path(self.model_dir, self.config_file)))

        # Remaining assertions are made here - each for atom data file, abundances
        # profile file and densities profile file.
        assert os.path.exists('/tmp/kurucz_cd23_chianti_H_He.h5'), \
            ('%s atom data file does not exist' % data_path(self.model_dir, self.config_file))

        assert os.path.exists(data_path('w7_13d', 'tardis_w7_13d_abundances.dat')), \
            ('%s abundances profile file does not exist' % data_path(self.model_dir, self.config_file))

        assert os.path.exists(data_path('w7_13d', 'tardis_w7_13d_densities.dat')), \
            ('%s densities profile file does not exist' % data_path(self.model_dir, self.model_dir))

        config_yaml['atom_data'] = '/tmp/kurucz_cd23_chianti_H_He.h5'

        if self.abundances is not None:
            config_yaml['model']['abundances']['filename'] = data_path(self.model_dir, self.abundances)

        if self.densities is not None:
            config_yaml['model']['structure']['filename'] = data_path(self.model_dir, self.densities)

        # The configuration hence obtained will be having appropriate file paths.
        tardis_config = Configuration.from_config_dict(config_yaml)

        # We now do a run with the prepared configuration and get radial1d model.
        self.obtained_w7_radial1d_model = Radial1DModel(tardis_config)
        simulation = Simulation(tardis_config)
        simulation.legacy_run_simulation(self.obtained_w7_radial1d_model)

        # The benchmark data against which assertions are to be made is ingested
        # from an already available compressed binary (.npz). This will return a
        # dictionary of numpy.ndarrays and the latter will also return dictionary
        # of numpy.ndarrays - difference is, they were astropy quantities earlier.
        self.expected_w7_ndarrays = np.load(os.path.join(data_path(self.model_dir, 'expected_ndarrays.npz')))
        self.expected_w7_astropy_quantities = np.load(os.path.join(data_path(self.model_dir, 'expected_quantities.npz')))

    def test_j_estimators(self):
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['j_estimators'],
            self.obtained_w7_radial1d_model.j_estimators)

    def test_j_blue_estimators(self):
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['j_blue_estimators'],
            self.obtained_w7_radial1d_model.j_blue_estimators)

        j_blues_norm_factor = self.expected_w7_astropy_quantities['j_blues_norm_factor']
        j_blues_norm_factor = j_blues_norm_factor * u.Unit('1 / (cm2 s)')
        np.testing.assert_allclose(j_blues_norm_factor,
                                   self.obtained_w7_radial1d_model.j_blues_norm_factor)

    def test_last_line_interactions(self):
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['last_line_interaction_in_id'],
            self.obtained_w7_radial1d_model.last_line_interaction_in_id)
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['last_line_interaction_out_id'],
            self.obtained_w7_radial1d_model.last_line_interaction_out_id)
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['last_line_interaction_shell_id'],
            self.obtained_w7_radial1d_model.last_line_interaction_shell_id)

        last_line_interaction_angstrom = self.expected_w7_astropy_quantities['last_line_interaction_angstrom']
        last_line_interaction_angstrom = last_line_interaction_angstrom * u.Unit('Angstrom')
        np.testing.assert_allclose(last_line_interaction_angstrom,
                                   self.obtained_w7_radial1d_model.last_line_interaction_angstrom)

    def test_nubar_estimators(self):
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['nubar_estimators'],
            self.obtained_w7_radial1d_model.nubar_estimators)

    def test_ws(self):
        np.testing.assert_allclose(
            self.expected_w7_ndarrays['ws'],
            self.obtained_w7_radial1d_model.ws)

    def test_spectrum(self):
        luminosity_inner = self.expected_w7_astropy_quantities['luminosity_inner']
        luminosity_inner = luminosity_inner * u.Unit('erg / s')
        np.testing.assert_allclose(luminosity_inner,
                                   self.obtained_w7_radial1d_model.luminosity_inner)

    def test_montecarlo_properties(self):
        montecarlo_luminosity = self.expected_w7_astropy_quantities['montecarlo_luminosity']
        montecarlo_luminosity = montecarlo_luminosity * u.Unit('erg / s')

        montecarlo_virtual_luminosity = self.expected_w7_astropy_quantities['montecarlo_virtual_luminosity']
        montecarlo_virtual_luminosity = montecarlo_virtual_luminosity * u.Unit('erg / s')

        montecarlo_nu = self.expected_w7_astropy_quantities['montecarlo_nu']
        montecarlo_nu = montecarlo_nu * u.Unit('Hz')

        np.testing.assert_allclose(montecarlo_luminosity,
                                   self.obtained_w7_radial1d_model.montecarlo_luminosity)

        np.testing.assert_allclose(montecarlo_virtual_luminosity,
                                   self.obtained_w7_radial1d_model.montecarlo_virtual_luminosity)

        np.testing.assert_allclose(montecarlo_nu,
                                   self.obtained_w7_radial1d_model.montecarlo_nu)

    def test_shell_temperature(self):
        t_rads = self.expected_w7_astropy_quantities['t_rads']
        t_rads = t_rads * u.Unit('K')

        np.testing.assert_allclose(t_rads, self.obtained_w7_radial1d_model.t_rads)
