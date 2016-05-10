import os
import yaml
import numpy as np
import pytest
from numpy.testing import assert_allclose
from astropy import units as u

from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


def data_path(fname):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "w7", fname)


@pytest.mark.skipif(not pytest.config.getoption("--slow"),
                    reason="slow tests can only be run using --slow")
@pytest.mark.skipif(not pytest.config.getvalue("baseline-data"),
                    reason="--baseline-data was not specified")
@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason="--atomic-dataset was not specified")
class TestW7:
    """
    Slow integration test for Stratified W7 setup.

    Assumed two compressed binaries (.npz) are placed in `baseline/w7`
    directory, whose path is provided by command line argument:

    * ndarrays.npz                       | * quantities.npz
    Contents (all (.npy)):               | Contents (all (.npy)):
        * last_interaction_type          |     * t_rads
        * last_line_interaction_out_id   |     * luminosity_inner
        * last_line_interaction_in_id    |     * montecarlo_luminosity
        * j_estimators                   |     * montecarlo_virtual_luminousity
        * j_blue_estimators              |     * time_of_simulation
        * last_line_interaction_shell_id |     * montecarlo_nu
        * nubar_estimators               |     * last_line_interaction_angstrom
        * ws                             |     * j_blues_norm_factor
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        This method does initial setup of creating configuration and performing
        a single run of integration test.
        """
        self.config_file = data_path("config_w7.yml")
        self.abundances = data_path("abundancies_w7.dat")
        self.densities = data_path("densities_w7.dat")

        # First we check whether atom data file exists at desired path.
        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
                                    pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(self.atom_data_filename), \
            "%s atom data file does not exist" % self.atom_data_filename

        # The available config file doesn't have file paths of atom data file,
        # densities and abundances profile files as desired. We form dictionary
        # from the config file and override those parameters by putting file
        # paths of these three files at proper places.
        config_yaml = yaml.load(open(self.config_file))
        config_yaml['atom_data'] = self.atom_data_filename
        config_yaml['model']['abundances']['filename'] = self.abundances
        config_yaml['model']['structure']['filename'] = self.densities

        # The config hence obtained will be having appropriate file paths.
        tardis_config = Configuration.from_config_dict(config_yaml)

        # We now do a run with prepared config and get radial1d model.
        self.obtained_radial1d_model = Radial1DModel(tardis_config)
        simulation = Simulation(tardis_config)
        simulation.legacy_run_simulation(self.obtained_radial1d_model)

        # The baseline data against which assertions are to be made is ingested
        # from already available compressed binaries (.npz). These will return
        # dictionaries of numpy.ndarrays for performing assertions.
        self.baseline_data_dir = os.path.join(os.path.expanduser(
                os.path.expandvars(pytest.config.getvalue('baseline-data'))), "w7")

        self.expected_ndarrays = np.load(os.path.join(self.baseline_data_dir,
                                                      "ndarrays.npz"))
        self.expected_quantities = np.load(os.path.join(self.baseline_data_dir,
                                                        "quantities.npz"))

    def test_j_estimators(self):
        assert_allclose(
                self.expected_ndarrays['j_estimators'],
                self.obtained_radial1d_model.j_estimators)

    def test_j_blue_estimators(self):
        assert_allclose(
                self.expected_ndarrays['j_blue_estimators'],
                self.obtained_radial1d_model.j_blue_estimators)

        j_blues_norm_factor = self.expected_quantities['j_blues_norm_factor']
        j_blues_norm_factor = j_blues_norm_factor * u.Unit('1 / (cm2 s)')

        assert_allclose(
                j_blues_norm_factor,
                self.obtained_radial1d_model.j_blues_norm_factor)

    def test_last_line_interactions(self):
        assert_allclose(
                self.expected_ndarrays['last_line_interaction_in_id'],
                self.obtained_radial1d_model.last_line_interaction_in_id)

        assert_allclose(
                self.expected_ndarrays['last_line_interaction_out_id'],
                self.obtained_radial1d_model.last_line_interaction_out_id)

        assert_allclose(
                self.expected_ndarrays['last_line_interaction_shell_id'],
                self.obtained_radial1d_model.last_line_interaction_shell_id)

        last_line_interaction_angstrom = self.expected_quantities['last_line_interaction_angstrom']
        last_line_interaction_angstrom = last_line_interaction_angstrom * u.Unit('Angstrom')

        assert_allclose(
                last_line_interaction_angstrom,
                self.obtained_radial1d_model.last_line_interaction_angstrom)

    def test_nubar_estimators(self):
        assert_allclose(
                self.expected_ndarrays['nubar_estimators'],
                self.obtained_radial1d_model.nubar_estimators)

    def test_ws(self):
        assert_allclose(
                self.expected_ndarrays['ws'],
                self.obtained_radial1d_model.ws)

    def test_spectrum(self):
        luminosity_inner = self.expected_quantities['luminosity_inner']
        luminosity_inner = luminosity_inner * u.Unit('erg / s')

        assert_allclose(
                luminosity_inner,
                self.obtained_radial1d_model.luminosity_inner)

    def test_montecarlo_properties(self):
        montecarlo_luminosity = self.expected_quantities['montecarlo_luminosity']
        montecarlo_luminosity = montecarlo_luminosity * u.Unit('erg / s')

        montecarlo_virtual_luminosity = self.expected_quantities['montecarlo_virtual_luminosity']
        montecarlo_virtual_luminosity = montecarlo_virtual_luminosity * u.Unit('erg / s')

        montecarlo_nu = self.expected_quantities['montecarlo_nu']
        montecarlo_nu = montecarlo_nu * u.Unit('Hz')

        assert_allclose(
                montecarlo_luminosity,
                self.obtained_radial1d_model.montecarlo_luminosity)

        assert_allclose(
                montecarlo_virtual_luminosity,
                self.obtained_radial1d_model.montecarlo_virtual_luminosity)

        assert_allclose(montecarlo_nu, self.obtained_radial1d_model.montecarlo_nu)

    def test_shell_temperature(self):
        t_rads = self.expected_quantities['t_rads']
        t_rads = t_rads * u.Unit('K')
        assert_allclose(t_rads, self.obtained_radial1d_model.t_rads)
