# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *
from astropy.tests.pytest_plugins import (
        pytest_addoption as _pytest_add_option
    )

import tardis
import pytest
from tardis.atomic import AtomData
from tardis.io.config_reader import Configuration
from tardis.io.util import yaml_load_config_file
from tardis.simulation import Simulation
from copy import deepcopy

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
    parser.addoption("--tardis-refdata", dest='tardis-refdata', default=None,
                     help="Path to Tardis Reference Folder")
    parser.addoption("--integration-tests",
                     dest="integration-tests", default=None,
                     help="path to configuration file for integration tests")
    parser.addoption("--generate-reference",
                     action="store_true", default=False,
                     help="execute integration test run to generate reference data")
    parser.addoption("--less-packets",
                     action="store_true", default=False,
                     help="Run integration tests with less packets.")

# -------------------------------------------------------------------------
# project specific fixtures
# -------------------------------------------------------------------------


@pytest.fixture(scope="session")
def atomic_data_fname():
    atomic_data_fname = os.path.join(pytest.config.getvalue(
        "tardis-refdata"), 'atom_data', 'kurucz_cd23_chianti_H_He.h5')
    if atomic_data_fname is None:
        pytest.skip('--tardis-refdata was not specified')
    else:
        return os.path.expandvars(os.path.expanduser(atomic_data_fname))

@pytest.fixture(scope="session")
def atomic_dataset(atomic_data_fname):
    atomic_data = AtomData.from_hdf5(atomic_data_fname)
    if atomic_data.md5 != '21095dd25faa1683f4c90c911a00c3f8':
        pytest.skip('Need default Kurucz atomic dataset '
                    '(md5="21095dd25faa1683f4c90c911a00c3f8"')
    else:
        return atomic_data

@pytest.fixture
def kurucz_atomic_data(atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
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

###
# HDF Fixtures
###

@pytest.fixture(scope="session")
def hdf_file_path(tmpdir_factory):
    path = tmpdir_factory.mktemp('hdf_buffer').join('test.hdf')
    return str(path)

@pytest.fixture(scope="session")
def config_verysimple():
    filename = 'tardis_configv1_verysimple.yml'
    path = os.path.abspath(os.path.join('tardis/io/tests/data/', filename))
    config = Configuration.from_yaml(path)
    return config

@pytest.fixture(scope="session")
def simulation_verysimple(config_verysimple, atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(4000)
    return sim
