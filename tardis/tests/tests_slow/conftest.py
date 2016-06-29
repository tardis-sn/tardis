import glob
import os
import yaml
import numpy as np
import pytest
from astropy import units as u

from tardis import __githash__ as tardis_githash
from tardis.tests.tests_slow.report import DokuReport
from tardis.tests.tests_slow.plot_helpers import PlotUploader


def pytest_addoption(parser):
    parser.addoption("--integration-tests", dest="integration-tests",
                     help="path to configuration file for integration tests")
    parser.addoption("--generate-reference", action="store_true",
                     help="execute integration test run to generate reference data")


def pytest_configure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath is not None:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        config.option.integration_tests_config = yaml.load(
            open(integration_tests_configpath))

        if not config.getoption("--generate-reference"):
            config.option.generate_reference['generate_reference'] = None
            # Used by DokuReport class to show build environment details in report.
            config._environment = []
            # prevent opening dokupath on slave nodes (xdist)
            if not hasattr(config, 'slaveinput'):
                config.dokureport = DokuReport(
                    config.option.integration_tests_config['dokuwiki'])
                config.pluginmanager.register(config.dokureport)


def pytest_unconfigure(config):
    # Unregister only if it was registered in pytest_configure
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath and not config.getoption("--generate-reference"):
        config.pluginmanager.unregister(config.dokureport)


@pytest.mark.hookwrapper
def pytest_runtest_makereport(item, call):
    # execute all other hooks to obtain the report object
    outcome = yield
    report = outcome.get_result()
    if report.when == "call":
        if "plot_object" in item.fixturenames:
            plot_obj = item.funcargs["plot_object"]
            plot_obj.upload(report)
            report.extra = plot_obj.get_extras()


@pytest.fixture(scope="function")
def plot_object(request):
    return PlotUploader(request)


@pytest.fixture(scope="class", params=[
    path for path in glob.glob(os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "*")) if os.path.isdir(path)
])
def data_path(request):
    integration_tests_config = request.config.option.integration_tests_config
    setup_name = os.path.basename(request.param)
    return {
        'config_dirpath': request.param,
        'reference_dirpath': os.path.join(os.path.expandvars(
            os.path.expanduser(integration_tests_config["reference"])), setup_name
        ),
        'setup_name': setup_name
    }


@pytest.fixture(scope="session")
def gen_ref_dirpath(request):
    integration_tests_config = request.config.integration_tests_config
    if request.config.getoption("--generate-reference"):
        path = os.path.join(os.path.expandvars(os.path.expanduser(
            integration_tests_config["generate_reference"])), tardis_githash
        )
        os.makedirs(os.path.join(path))
    else:
        path = None
    return path


@pytest.fixture(scope="class")
def reference(request, data_path):
    """
    Fixture to ingest reference data for slow test from already available
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
    reference_dirpath = data_path['reference_dirpath']

    # TODO: make this fixture ingest data from an HDF5 file.
    ndarrays = dict(np.load(os.path.join(reference_dirpath, "ndarrays.npz")))
    quantities = dict(np.load(os.path.join(reference_dirpath, "quantities.npz")))
    spectrum = dict(np.load(os.path.join(reference_dirpath, "spectrum.npz")))

    # Associate CGS units to ndarrays of reference quantities.
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
