import os
import shutil
import tempfile
import yaml
import numpy as np
import pytest
from astropy import units as u

from tardis.tests.tests_slow.report import DokuReport


def pytest_configure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath is not None:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        config.option.integration_tests_config = yaml.load(
            open(integration_tests_configpath))

        # Used by DokuReport class to show build environment details in report.
        config._environment = []
        # prevent opening dokupath on slave nodes (xdist)
        if not hasattr(config, 'slaveinput'):
            config.dokureport = DokuReport(
                config.option.integration_tests_config['dokuwiki'])
            config.pluginmanager.register(config.dokureport)
    # A common tempdir for storing plots / PDFs and other slow test related data
    # generated during execution.
    tempdir_session = tempfile.mkdtemp()
    config.option.tempdir = tempdir_session


def pytest_unconfigure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath is not None:
        config.pluginmanager.unregister(config.dokureport)
    # Remove tempdir by recursive deletion
    shutil.rmtree(config.option.tempdir)


@pytest.fixture(scope="session")
def integration_tests_config(request):
    return request.config.option.integration_tests_config


@pytest.fixture(scope="session")
def reference_datadir(integration_tests_config):
    return os.path.expandvars(
        os.path.expanduser(integration_tests_config["reference"])
    )


@pytest.fixture(scope="session")
def data_path():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "w7")


@pytest.fixture(scope="session")
def reference(request, reference_datadir):
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

    # TODO: make this fixture ingest data from an HDF5 file.
    ndarrays = dict(np.load(os.path.join(reference_datadir, "ndarrays.npz")))
    quantities = dict(np.load(os.path.join(reference_datadir, "quantities.npz")))
    spectrum = dict(np.load(os.path.join(reference_datadir, "spectrum.npz")))

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
