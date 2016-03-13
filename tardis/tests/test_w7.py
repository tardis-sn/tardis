import os
import yaml
import pytest
import numpy as np
import numpy.testing as nptesting
from astropy import units as u

import tardis
from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


def w7_path(fname):
    return os.path.join(tardis.__path__[0], 'tests', 'w7_13d', fname)


@pytest.mark.skipif(not pytest.config.getoption("--run-slow"),
                    reason='this is a slow test, add --run-slow to run')
class TestW7:
    """
    Integration test for a run with the Stratified W7 setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        This method does an initial setup of making the configuration and performing
        a single run of this integration test from Stratified W7 setup.
        """

        # First we check whether the config file does even exist at the desired path.
        assert os.path.exists(w7_path('tardis_w7.yml')), ('%s config file does not exist' %
                                                          w7_path('tardis_w7.yml'))

        # The available config file does not have the filepaths of atom data file,
        # densities and abundances profile files as desired. We form a dictionary
        # from the config file and tweak those parameters - putting the filepaths
        # of these three files at proper places.
        self.config_yaml = yaml.load(open(w7_path('tardis_w7.yml')))

        # Remaining assertions are made here - each for atom data file, abundances
        # profile file and densities profile file.
        assert os.path.exists('/tmp/kurucz_cd23_chianti_H_He.h5'), ('%s atom data file does not exist' %
                                                                    w7_path('tardis_w7.yml'))

        assert os.path.exists(w7_path('tardis_w7_13d_abundances.dat')), ('%s abundances profile file does not exist' %
                                                                         w7_path('tardis_w7.yml'))

        assert os.path.exists(w7_path('tardis_w7_13d_densities.dat')), ('%s densities profile file does not exist' %
                                                                        w7_path('tardis_w7.yml'))

        self.config_yaml['atom_data'] = '/tmp/kurucz_cd23_chianti_H_He.h5'
        self.config_yaml['model']['abundances']['filename'] = w7_path('tardis_w7_13d_abundances.dat')
        self.config_yaml['model']['structure']['filename'] = w7_path('tardis_w7_13d_densities.dat')

        # The configuration hence obtained will be having appropriate file paths.
        tardis_config = Configuration.from_config_dict(self.config_yaml)

        # We now do a run with the prepared configuration and get radial1d model.
        self.obtained_w7_radial1d_model = Radial1DModel(tardis_config)
        self.simulation = Simulation(tardis_config)
        self.simulation.legacy_run_simulation(self.obtained_w7_radial1d_model)

        # The benchmark data against which assertions are to be made is ingested
        # from an already available compressed binary (.npz). This will return a
        # dictionary of numpy.ndarrays and the latter will also return dictionary
        # of numpy.ndarrays - difference is, they were astropy quantities earlier.
        self.expected_w7_ndarrays = np.load(os.path.join(w7_path('expected_ndarrays.npz')))
        self.expected_w7_astropy_quantities = np.load(os.path.join(w7_path('expected_quantities.npz')))

    def test_j_estimators(self):
        np.testing.assert_allclose(
                self.expected_w7_ndarrays['j_estimators'],
                self.obtained_w7_radial1d_model.j_estimators)

    def test_j_blue_estimators(self):
        np.testing.assert_allclose(
                self.expected_w7_ndarrays['j_blue_estimators'],
                self.obtained_w7_radial1d_model.j_blue_estimators)

        self.expected_w7_astropy_quantities['j_blues_norm_factor'] *= u.Unit('1 / (cm2 s)')
        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['j_blues_norm_factor'],
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

        self.expected_w7_astropy_quantities['last_line_interaction_angstrom'] *= u.Unit('Angstrom')
        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['last_line_interaction_angstrom'],
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
        self.expected_w7_astropy_quantities['luminosity_inner'] *= u.Unit('erg / s')
        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['luminosity_inner'],
                self.obtained_w7_radial1d_model.luminosity_inner)

    def test_montecarlo_properties(self):
        self.expected_w7_astropy_quantities['montecarlo_luminosity'] *= u.Unit('erg / s')
        self.expected_w7_astropy_quantities['montecarlo_virtual_luminosity'] *= u.Unit('erg / s')
        self.expected_w7_astropy_quantities['montecarlo_nu'] *= u.Unit('Hz')

        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['montecarlo_luminosity'],
                self.obtained_w7_radial1d_model.montecarlo_luminosity)

        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['montecarlo_virtual_luminosity'],
                self.obtained_w7_radial1d_model.montecarlo_virtual_luminosity)

        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['montecarlo_nu'],
                self.obtained_w7_radial1d_model.montecarlo_nu)

    def test_shell_temperature(self):
        self.expected_w7_astropy_quantities['t_rads'] *= u.Unit('K')

        np.testing.assert_allclose(
                self.expected_w7_astropy_quantities['t_rads'],
                self.obtained_w7_radial1d_model.t_rads)
