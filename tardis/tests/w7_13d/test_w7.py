import pytest
import numpy as np
import tardis
import numpy.testing as nptesting
from astropy import units as u
import os


@pytest.mark.skipif(not pytest.config.getoption("--run-slow"),
                    reason='this is a slow test, add --run-slow to run')
class TestW7:
    """
    Integration test for a run with the Stratified W7 setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        self.w7_config_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'tardis_w7.yml')
        self.atom_data_filepath = os.path.realpath('/tmp/kurucz_cd23_chianti_H_He.h5')
        self.obtained_w7_radial1d_model = tardis.run_tardis(config=self.w7_config_filepath,
                                                            atom_data=self.atom_data_filepath)

        self.expected_w7_ndarrays = np.load(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 'expected_ndarrays.npz')
        )

        self.expected_w7_astropy_quantities = np.load(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 'expected_quantities.npz')
        )

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
