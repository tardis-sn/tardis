# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.
import copy
import shutil
import tempfile
import zipfile

from astropy.tests.helper import remote_data
from astropy.tests.pytest_plugins import *
from astropy.tests.pytest_plugins import (
        pytest_addoption as _pytest_add_option
    )
from astropy.utils.data import download_file

import tardis
from tardis.atomic import AtomData

###
# Astropy
###

# Uncomment the following line to treat all DeprecationWarnings as
# exceptions
# enable_deprecations_as_exceptions()

# Uncomment and customize the following lines to add/remove entries from
# the list of packages for which version numbers are displayed when running
# the tests. Making it pass for KeyError is essential in some cases when
# the package uses other astropy affiliated packages.
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

# Uncomment the following lines to display the version number of the
# package rather than the version number of Astropy in the top line when
# running the tests.
import os

# This is to figure out the affiliated package version, rather than
# using Astropy's
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
    parser.addoption("--integration-tests", dest="integration-tests", default=None,
                     help="path to configuration file for integration tests")
    parser.addoption("--generate-reference", action="store_true", default=False,
                     help="execute integration test run to generate reference data")
    parser.addoption("--less-packets", action="store_true", default=False,
                     help="Run integration tests with less packets.")


# -------------------------------------------------------------------------
# project specific fixtures
# -------------------------------------------------------------------------


@pytest.fixture
def test_data_path():
    return os.path.join(tardis.__path__[0], 'tests', 'data')


@pytest.fixture
def included_he_atomic_data(test_data_path):
    atomic_db_fname = os.path.join(test_data_path, 'chianti_he_db.h5')
    return AtomData.from_hdf5(atomic_db_fname)


@pytest.fixture(scope="session")
def tardis_config_verysimple():
    return 'tardis/io/tests/data/tardis_configv1_verysimple.yml'


@pytest.fixture(scope="session")
def atom_data(request):
    atom_data_name = 'kurucz_cd23_chianti_H_He.h5'
    # Download and cache the zip file of atomic data.
    atom_data_cache = download_file(
        'http://www.mpa-garching.mpg.de/~michi/tardis/data/kurucz_cd23_chianti_H_He.zip',
        cache=True
    )

    # Obtained file is a zip file, hence unzipped inside a tempdir.
    atom_data_zipfile = zipfile.ZipFile(atom_data_cache)
    atom_data_extract_tempdir = tempfile.mkdtemp()
    atom_data_zipfile.extract(atom_data_name, path=atom_data_extract_tempdir)

    atom_data = AtomData.from_hdf5(
        os.path.join(atom_data_extract_tempdir, atom_data_name)
    )

    if atom_data.md5 != '21095dd25faa1683f4c90c911a00c3f8':
        pytest.skip('Need default Kurucz atomic dataset '
                    '(md5="21095dd25faa1683f4c90c911a00c3f8"')

    def fin():
        # Delete the tempdir as no longer required.
        shutil.rmtree(atom_data_extract_tempdir)
    request.addfinalizer(fin)

    return atom_data
