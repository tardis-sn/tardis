# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *
from astropy.tests.pytest_plugins import (
        pytest_addoption as _pytest_add_option,
        pytest_configure as _pytest_configure,
        pytest_unconfigure as _pytest_unconfigure
        )

import yaml
import tempfile
import tardis
import pytest
from tardis.atomic import AtomData
from tardis.io.config_reader import Configuration
from tardis.io.util import yaml_load_config_file

# For specifying error while exception handling
from socket import gaierror

try:
    import dokuwiki
except ImportError:
    dokuwiki_available = False
else:
    dokuwiki_available = True


###
# Astropy
###

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions
# enable_deprecations_as_exceptions()

## Uncomment and customize the following lines to add/remove entries from
## the list of packages for which version numbers are displayed when running
## the tests. Making it pass for KeyError is essential in some cases when
## the package uses other astropy affiliated packages.
try:
    PYTEST_HEADER_MODULES['Numpy'] = 'numpy'
    PYTEST_HEADER_MODULES['Scipy'] = 'scipy'
    PYTEST_HEADER_MODULES['Pandas'] = 'pandas'
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['Yaml'] = 'yaml'
    PYTEST_HEADER_MODULES['Cython'] = 'cython'
    PYTEST_HEADER_MODULES['h5py'] = 'h5py'
    PYTEST_HEADER_MODULES['Matplotlib'] = 'matplotlib'
    PYTEST_HEADER_MODULES['Ipython'] = 'IPython'
#     del PYTEST_HEADER_MODULES['h5py']
except (NameError, KeyError):  # NameError is needed to support Astropy < 1.0
    pass

## Uncomment the following lines to display the version number of the
## package rather than the version number of Astropy in the top line when
## running the tests.
import os

## This is to figure out the affiliated package version, rather than
## using Astropy's
try:
    from .version import version
except ImportError:
    version = 'dev'

try:
    packagename = os.path.basename(os.path.dirname(__file__))
    TESTED_VERSIONS[packagename] = version
except NameError:   # Needed to support Astropy <= 1.0.0
    pass

# -------------------------------------------------------------------------
# Initialization
# -------------------------------------------------------------------------


def pytest_addoption(parser):
    _pytest_add_option(parser)
    parser.addoption("--atomic-dataset", dest='atomic-dataset', default=None,
                     help="filename for atomic dataset")

    parser.addoption("--integration-tests", dest="integration-tests",
                     help="path to configuration file for integration tests")


def pytest_configure(config):
    _pytest_configure(config)
    html_file = tempfile.NamedTemporaryFile(delete=False)
    # Html test report will be generated at this filepath by pytest-html plugin
    config.option.htmlpath = html_file.name


def pytest_unconfigure(config):
    _pytest_unconfigure()
    # Html report created by pytest-html plugin is read here, uploaded to
    # dokuwiki and finally deleted.
    if dokuwiki_available:
        githash = tardis.__githash__
        report_content = open(config.option.htmlpath, 'rb').read()
        report_content = report_content.replace("<!DOCTYPE html>", "")

        report_content = (
            "Test executed on commit "
            "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n"
            "{1}".format(githash, report_content)
        )

        try:
            doku_conn = dokuwiki.DokuWiki(
                url=integration_tests_config["dokuwiki"]["url"],
                user=integration_tests_config["dokuwiki"]["username"],
                password=integration_tests_config["dokuwiki"]["password"])
        except gaierror, dokuwiki.DokuWikiError:
            print "Dokuwiki connection not established, report upload failed!"
        else:
            # Upload report on dokuwiki. Temporary link due to prototyping purposes.
            doku_conn.pages.append("reports:{0}".format(githash[:7]), report_content)
            print "Uploaded report on Dokuwiki."

    # Remove the local report file. Keeping the report saved on local filesystem
    # is not desired, hence deleted.
    os.unlink(config.option.htmlpath)
    print "Deleted temporary file containing html report."

# -------------------------------------------------------------------------
# hooks for influencing reporting (invoked from _pytest_terminal)
# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
# project specific fixtures
# -------------------------------------------------------------------------


@pytest.fixture(scope="session")
def atomic_data_fname():
    atomic_data_fname = pytest.config.getvalue("atomic-dataset")
    if atomic_data_fname is None:
        pytest.skip('--atomic_database was not specified')
    else:
        return os.path.expandvars(os.path.expanduser(atomic_data_fname))


@pytest.fixture
def kurucz_atomic_data(atomic_data_fname):
    atomic_data = AtomData.from_hdf5(atomic_data_fname)

    if atomic_data.md5 != '21095dd25faa1683f4c90c911a00c3f8':
        pytest.skip('Need default Kurucz atomic dataset '
                    '(md5="21095dd25faa1683f4c90c911a00c3f8"')
    else:
        return atomic_data

@pytest.fixture
def test_data_path():
    return os.path.join(tardis.__path__[0], 'tests', 'data')

@pytest.fixture
def included_he_atomic_data(test_data_path):
    atomic_db_fname = os.path.join(test_data_path, 'chianti_he_db.h5')
    return AtomData.from_hdf5(atomic_db_fname)

@pytest.fixture
def tardis_config_verysimple():
    return yaml_load_config_file(
        'tardis/io/tests/data/tardis_configv1_verysimple.yml')


@pytest.fixture(scope="session")
def integration_tests_config(request):
    integration_tests_configpath = request.config.getvalue("integration-tests")
    if integration_tests_configpath is None:
        pytest.skip('--integration-tests was not specified')
    else:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        return yaml.load(open(integration_tests_configpath))
