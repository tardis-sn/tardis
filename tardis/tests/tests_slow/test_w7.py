import os
import yaml
import numpy as np
import pytest
from numpy.testing import assert_allclose
from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u

from tardis.atomic import AtomData
from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration
from tardis.tests.tests_slow import runslow

@pytest.fixture(scope="module")
def data_path():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "w7")


@pytest.fixture(scope="module")
def baseline(request):
    """
    Fixture to ingest baseline data for slow test from already available
    compressed binaries (.npz). All data is collected in one dict and
    returned away.

    Assumed three compressed binaries (.npz) are placed in `slow-test-data/w7`
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

    * spectrum.npz
    Contents (all (.npy)):
        * luminosity_density_nu          | * wavelength
        * delta_frequency                | * luminosity_density_lambda
    """

    # TODO: make this fixture ingest data from an HDF5 file.
    datadir_path = os.path.join(os.path.expanduser(os.path.expandvars(
            request.config.getvalue('slow-test-data'))), "w7")

    ndarrays = dict(np.load(os.path.join(datadir_path, "ndarrays.npz")))
    quantities = dict(np.load(os.path.join(datadir_path, "quantities.npz")))
    spectrum = dict(np.load(os.path.join(datadir_path, "spectrum.npz")))

    # Associate CGS units to ndarrays of baseline quantities.
    ndarrays.update(
        j_blues_norm_factor=
            u.Quantity(quantities['j_blues_norm_factor'], '1 / (cm2 s)'),
        last_line_interaction_angstrom=
            u.Quantity(quantities['last_line_interaction_angstrom'], 'Angstrom'),
        luminosity_inner=
            u.Quantity(quantities['luminosity_inner'], 'erg / s'),
        montecarlo_luminosity=
            u.Quantity(quantities['montecarlo_luminosity'], 'erg / s'),
        montecarlo_virtual_luminosity=
            u.Quantity(quantities['montecarlo_virtual_luminosity'], 'erg / s'),
        montecarlo_nu=
            u.Quantity(quantities['montecarlo_nu'], 'Hz'),
        t_rads=
            u.Quantity(quantities['t_rads'], 'K'),
        luminosity_density_nu=
            u.Quantity(spectrum['luminosity_density_nu'], 'erg'),
        delta_frequency=
            u.Quantity(spectrum['delta_frequency'], 'Hz'),
        wavelength=
            u.Quantity(spectrum['wavelength'], 'Angstrom'),
        luminosity_density_lambda=
            u.Quantity(spectrum['luminosity_density_lambda'], 'erg / (Angstrom s)'))
    return ndarrays


@runslow
class TestW7(object):
    """
    Slow integration test for Stratified W7 setup.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, baseline, data_path):
        """
        This method does initial setup of creating configuration and performing
        a single run of integration test.
        """
        self.config_file = os.path.join(data_path, "config_w7.yml")
        self.abundances = os.path.join(data_path, "abundancies_w7.dat")
        self.densities = os.path.join(data_path, "densities_w7.dat")

        # First we check whether atom data file exists at desired path.
        self.atom_data_filename = os.path.expanduser(os.path.expandvars(
                                    pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(self.atom_data_filename), \
            "{0} atom data file does not exist".format(self.atom_data_filename)

        # The available config file doesn't have file paths of atom data file,
        # densities and abundances profile files as desired. We load the atom
        # data seperately and provide it to tardis_config later. For rest of
        # the two, we form dictionary from the config file and override those
        # parameters by putting file paths of these two files at proper places.
        config_yaml = yaml.load(open(self.config_file))
        config_yaml['model']['abundances']['filename'] = self.abundances
        config_yaml['model']['structure']['filename'] = self.densities

        # Load atom data file separately, pass it for forming tardis config.
        self.atom_data = AtomData.from_hdf5(self.atom_data_filename)

        # Check whether the atom data file in current run and the atom data
        # file used in obtaining the baseline data for slow tests are same.
        # TODO: hard coded UUID for kurucz atom data file, generalize it later.
        kurucz_data_file_uuid1 = "5ca3035ca8b311e3bb684437e69d75d7"
        assert self.atom_data.uuid1 == kurucz_data_file_uuid1

        # The config hence obtained will be having appropriate file paths.
        tardis_config = Configuration.from_config_dict(config_yaml, self.atom_data)

        # We now do a run with prepared config and get radial1d model.
        self.obtained_radial1d_model = Radial1DModel(tardis_config)
        simulation = Simulation(tardis_config)
        simulation.legacy_run_simulation(self.obtained_radial1d_model)

        # Get the baseline data through the fixture.
        self.baseline = baseline

    def test_j_estimators(self):
        assert_allclose(
                self.baseline['j_estimators'],
                self.obtained_radial1d_model.j_estimators)

    def test_j_blue_estimators(self):
        assert_allclose(
                self.baseline['j_blue_estimators'],
                self.obtained_radial1d_model.j_blue_estimators)

        assert_quantity_allclose(
                self.baseline['j_blues_norm_factor'],
                self.obtained_radial1d_model.j_blues_norm_factor)

    def test_last_line_interactions(self):
        assert_allclose(
                self.baseline['last_line_interaction_in_id'],
                self.obtained_radial1d_model.last_line_interaction_in_id)

        assert_allclose(
                self.baseline['last_line_interaction_out_id'],
                self.obtained_radial1d_model.last_line_interaction_out_id)

        assert_allclose(
                self.baseline['last_line_interaction_shell_id'],
                self.obtained_radial1d_model.last_line_interaction_shell_id)

        assert_quantity_allclose(
                self.baseline['last_line_interaction_angstrom'],
                self.obtained_radial1d_model.last_line_interaction_angstrom)

    def test_nubar_estimators(self):
        assert_allclose(
                self.baseline['nubar_estimators'],
                self.obtained_radial1d_model.nubar_estimators)

    def test_ws(self):
        assert_allclose(
                self.baseline['ws'],
                self.obtained_radial1d_model.ws)

    def test_luminosity_inner(self):
        assert_quantity_allclose(
                self.baseline['luminosity_inner'],
                self.obtained_radial1d_model.luminosity_inner)

    def test_spectrum(self):
        assert_quantity_allclose(
                self.baseline['luminosity_density_nu'],
                self.obtained_radial1d_model.spectrum.luminosity_density_nu)

        assert_quantity_allclose(
                self.baseline['delta_frequency'],
                self.obtained_radial1d_model.spectrum.delta_frequency)

        assert_quantity_allclose(
                self.baseline['wavelength'],
                self.obtained_radial1d_model.spectrum.wavelength)

        assert_quantity_allclose(
                self.baseline['luminosity_density_lambda'],
                self.obtained_radial1d_model.spectrum.luminosity_density_lambda)

    def test_montecarlo_properties(self):
        assert_quantity_allclose(
                self.baseline['montecarlo_luminosity'],
                self.obtained_radial1d_model.montecarlo_luminosity)

        assert_quantity_allclose(
                self.baseline['montecarlo_virtual_luminosity'],
                self.obtained_radial1d_model.montecarlo_virtual_luminosity)

        assert_quantity_allclose(
                self.baseline['montecarlo_nu'],
                self.obtained_radial1d_model.montecarlo_nu)

    def test_shell_temperature(self):
        assert_quantity_allclose(
                self.baseline['t_rads'],
                self.obtained_radial1d_model.t_rads)
