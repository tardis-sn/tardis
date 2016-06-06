import os
import shutil
import tempfile
import yaml
import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import remote_data
import tardis
from tardis.tests.tests_slow.report import DokuReport

# For specifying error while exception handling
from socket import gaierror

try:
    import dokuwiki
except ImportError:
    dokuwiki_available = False
else:
    dokuwiki_available = True


def pytest_configure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath is not None:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        config.option.integration_tests_config = yaml.load(
            open(integration_tests_configpath))

        dokufile = tempfile.NamedTemporaryFile(delete=False)
        # Report will be generated at this filepath by pytest-html plugin
        config.option.dokufile = dokufile

        # prevent opening dokupath on slave nodes (xdist)
        if not hasattr(config, 'slaveinput'):
            config.dokureport = DokuReport(config.option.dokufile)
            config.pluginmanager.register(config.dokureport)
    # A common tempdir for storing plots / PDFs and other slow test related data
    # generated during execution.
    tempdir_session = tempfile.mkdtemp()
    config.option.tempdir = tempdir_session


@remote_data
def pytest_unconfigure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath is not None:
        # Html report created by pytest-html plugin is read here, uploaded to
        # dokuwiki and finally deleted.
        if dokuwiki_available:
            githash = tardis.__githash__
            report_content = open(config.option.dokufile.name, 'rb').read()
            report_content = report_content.replace("<!DOCTYPE html>", "")

            report_content = (
                "Test executed on commit "
                "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n"
                "{1}".format(githash, report_content)
            )

            # These steps are already performed by `integration_tests_config` but
            # all the fixtures are teared down and no longer usable, when this
            # method is being called by pytest, hence they are called as functions.
            integration_tests_config = config.option.integration_tests_config
            try:
                doku_conn = dokuwiki.DokuWiki(
                    url=integration_tests_config["dokuwiki"]["url"],
                    user=integration_tests_config["dokuwiki"]["username"],
                    password=integration_tests_config["dokuwiki"]["password"])
            except gaierror, dokuwiki.DokuWikiError:
                print "Dokuwiki connection not established, report upload failed!"
            else:
                # Upload report on dokuwiki. Temporary link due to prototyping purposes.
                doku_conn.pages.set("reports:{0}".format(githash[:7]), report_content)
                print "Uploaded report on Dokuwiki."
        config.pluginmanager.unregister(config.dokureport)

    # Remove the local report file. Keeping the report saved on local filesystem
    # is not desired, hence deleted.
    os.unlink(config.option.dokufile.name)
    print "Deleted temporary file containing html report."
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
