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
from tardis.montecarlo.base import MontecarloRunner
from tardis.plasma.standard_plasmas import LegacyPlasmaArray

def data_path(fname):
    return os.path.join(tardis.__path__[0], 'tests', 'data', fname)


path_abn = "/tardis/io/tests/data/abn_tom_test.yml"
path_y7 = "/tardis/io/tests/data/tardis_w7.yml"
path_simple = "/tardis/io/tests/data/tardis_configv1_verysimple.yml"

slow = pytest.mark.skipif( not pytest.config.getoption('--slow') , reason= "needs --slow to run this")


@slow
@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
@pytest.fixture(scope="module", autouse=True, params=[path_abn, path_y7,path_simple])
def setup(request):
    """
    Parametrized fixture. Test will be run with all parameter specified in "params" array
    -------

    """
    atom_data_filename = os.path.expanduser(os.path.expandvars(
            pytest.config.getvalue('atomic-dataset')))
    assert os.path.exists(atom_data_filename), ("{0} atomic datafiles"
                                                         " does not seem to "
                                                         "exist".format(atom_data_filename))

    model = Initiate_slow_configs(request.param, atom_data_filename)
    model.simulation.legacy_run_simulation(setup.model)

    return model

def test_spectrum(setup):
    if setup.path == path_abn:
        luminosity_density = np.load(data_path("abn_tom_test_luminosity_density_lambda.npy"))
    elif setup.path == path_y7:
        luminosity_density = np.load(data_path("w7_luminosity_density_lambda.npy"))
    elif setup.path == path_simple:
         luminosity_density = np.load(data_path("simple_test_luminosity_density_lambda.npy"))

    luminosity_density = luminosity_density * u.Unit(
        'erg / (Angstrom s)')

    np.testing.assert_allclose(
        setup.model.spectrum.luminosity_density_lambda,luminosity_density)

def test_virtual_spectrum(setup):
    if setup.path == path_abn:
        virtual_luminosity_density = np.load(
            data_path("abn_tom_test_virtual_luminosity_density_lambda.npy"))
    elif setup.path == path_y7:
        virtual_luminosity_density = np.load(data_path("w7_virtual_luminosity_density_lambda.npy"))
    elif setup.path == path_simple:
         virtual_luminosity_density = np.load(data_path("simple_test_virtual_luminosity_density_lambda.npy"))
    virtual_luminosity_density = virtual_luminosity_density * u.Unit(
            'erg / (Angstrom s)')
    np.testing.assert_allclose(
            setup.model.spectrum_virtual.luminosity_density_lambda,
            virtual_luminosity_density)

def test_j_blue_estimators(setup):

    if setup.path == path_abn:
        j_blue_estimator = np.load(data_path("abn_tom_jblues_estimator.npy"))
    elif setup.path == path_y7:
        j_blue_estimator = np.load(data_path("w7_jblues_estimator.npy"))
    elif setup.path == path_simple:
        j_blue_estimator = np.load(data_path("simple_test_j_blue_estimator.npy"))
    if j_blue_estimator == None:
        raise Exception('J_blue file not found')
    np.testing.assert_allclose(setup.model.runner.j_blue_estimator, j_blue_estimator)
def test_runner_properties(setup):
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

        assert type(getattr(setup.model.runner, prop)) == prop_type, ("wrong type of attribute '{}': expected {}, found {}".format(prop, prop_type, type(getattr(setup.model.runner, prop))))


def test_legacy_model_properties(setup):
    """Tests whether a number of model attributes exist and also verifies
    their types

    Currently, model attributes needed to run the gui and to call the model
    routine to_hdf5 are checked.

    Notes
    -----
    The list of properties may be incomplete

    """

    props_required_by_gui = dict([
            ("converged", bool),
            ("iterations_executed", int),
            ("iterations_max_requested", int),
            ("current_no_of_packets", int),
            ("no_of_packets", int),
            ("no_of_virtual_packets", int),
            ])
    props_required_by_tohdf5 = dict([
            ("runner", MontecarloRunner),
            ("plasma_array", LegacyPlasmaArray),
            ("last_line_interaction_in_id", np.ndarray),
            ("last_line_interaction_out_id", np.ndarray),
            ("last_line_interaction_shell_id", np.ndarray),
            ("last_line_interaction_in_id", np.ndarray),
            ("last_line_interaction_angstrom", u.quantity.Quantity),
            ])

    required_props = props_required_by_gui.copy()
    required_props.update(props_required_by_tohdf5)

    for prop, prop_type in required_props.items():

        assert type(getattr(setup.model, prop)) == prop_type, ("wrong type of attribute '{}': expected {}, found {}".format(prop, prop_type, type(getattr(setup.model, prop))))



class Initiate_slow_configs():
    """
    Initiate class with tardis configs, that will be passed
    as an argument to test function
    """
    def __init__(self,path, atomic_dataset):
        self.path = path
        self.config_yaml = yaml.load(open(path))
        self.config_yaml['atom_data'] = atomic_dataset

        tardis_config = Configuration.from_config_dict(self.config_yaml)

        self.model = Radial1DModel(tardis_config)
        self.simulation = Simulation(tardis_config)

