# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.

from astropy.version import version as astropy_version

if astropy_version < "3.0":
    # With older versions of Astropy, we actually need to import the pytest
    # plugins themselves in order to make them discoverable by pytest.
    from astropy.tests.pytest_plugins import *
else:
    # As of Astropy 3.0, the pytest plugins provided by Astropy are
    # automatically made available when Astropy is installed. This means it's
    # not necessary to import them here, but we still need to import global
    # variables that are used for configuration.
    from astropy.tests.plugins.display import (
        PYTEST_HEADER_MODULES,
        TESTED_VERSIONS,
    )

from astropy.tests.helper import enable_deprecations_as_exceptions

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions. For Astropy v2.0 or later, there are 2 additional keywords,
## as follow (although default should work for most cases).
## To ignore some packages that produce deprecation warnings on import
## (in addition to 'compiler', 'scipy', 'pygments', 'ipykernel', and
## 'setuptools'), add:
##     modules_to_ignore_on_import=['module_1', 'module_2']
## To ignore some specific deprecation warning messages for Python version
## MAJOR.MINOR or later, add:
##     warnings_to_ignore_by_pyver={(MAJOR, MINOR): ['Message to ignore']}
# enable_deprecations_as_exceptions()

## Uncomment and customize the following lines to add/remove entries from
## the list of packages for which version numbers are displayed when running
## the tests. Making it pass for KeyError is essential in some cases when
## the package uses other astropy affiliated packages.
# try:
#     PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
#     PYTEST_HEADER_MODULES['scikit-image'] = 'skimage'
#     del PYTEST_HEADER_MODULES['h5py']
# except (NameError, KeyError):  # NameError is needed to support Astropy < 1.0
#     pass

## Uncomment the following lines to display the version number of the
## package rather than the version number of Astropy in the top line when
## running the tests.
# import os
#
## This is to figure out the package version, rather than
## using Astropy's
# try:
#     from .version import version
# except ImportError:
#     version = 'dev'
#
# try:
#     packagename = os.path.basename(os.path.dirname(__file__))
#     TESTED_VERSIONS[packagename] = version
# except NameError:   # Needed to support Astropy <= 1.0.0
#     pass

### XXX #### Here the TARDIS testing stuff begins ### XXX ####

import pytest

from tardis.io.util import yaml_load_config_file
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
import pandas as pd

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
    PYTEST_HEADER_MODULES["Numpy"] = "numpy"
    PYTEST_HEADER_MODULES["Scipy"] = "scipy"
    PYTEST_HEADER_MODULES["Pandas"] = "pandas"
    PYTEST_HEADER_MODULES["Astropy"] = "astropy"
    PYTEST_HEADER_MODULES["Yaml"] = "yaml"
    PYTEST_HEADER_MODULES["Cython"] = "cython"
    PYTEST_HEADER_MODULES["h5py"] = "h5py"
    PYTEST_HEADER_MODULES["Matplotlib"] = "matplotlib"
    PYTEST_HEADER_MODULES["Ipython"] = "IPython"
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
    version = "dev"

try:
    packagename = os.path.basename(os.path.dirname(__file__))
    TESTED_VERSIONS[packagename] = version
except NameError:  # Needed to support Astropy <= 1.0.0
    pass


# -------------------------------------------------------------------------
# Initialization
# -------------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--tardis-refdata", default=None, help="Path to Tardis Reference Folder"
    )
    parser.addoption(
        "--integration-tests",
        dest="integration-tests",
        default=None,
        help="path to configuration file for integration tests",
    )
    parser.addoption(
        "--generate-reference",
        action="store_true",
        default=False,
        help="generate reference data instead of testing",
    )
    parser.addoption(
        "--less-packets",
        action="store_true",
        default=False,
        help="Run integration tests with less packets.",
    )


# -------------------------------------------------------------------------
# project specific fixtures
# -------------------------------------------------------------------------


@pytest.fixture(scope="session")
def generate_reference(request):
    option = request.config.getoption("--generate-reference")
    if option is None:
        return False
    else:
        return option


@pytest.fixture(scope="session")
def tardis_ref_path(request):
    tardis_ref_path = request.config.getoption("--tardis-refdata")
    if tardis_ref_path is None:
        pytest.skip("--tardis-refdata was not specified")
    else:
        return os.path.expandvars(os.path.expanduser(tardis_ref_path))


from tardis.tests.fixtures.atom_data import *


@pytest.yield_fixture(scope="session")
def tardis_ref_data(tardis_ref_path, generate_reference):
    if generate_reference:
        mode = "w"
    else:
        mode = "r"
    with pd.HDFStore(
        os.path.join(tardis_ref_path, "unit_test_data.h5"), mode=mode
    ) as store:
        yield store


@pytest.fixture
def tardis_config_verysimple():
    return yaml_load_config_file(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )


###
# HDF Fixtures
###


@pytest.fixture(scope="session")
def hdf_file_path(tmpdir_factory):
    path = tmpdir_factory.mktemp("hdf_buffer").join("test.hdf")
    return str(path)


@pytest.fixture(scope="session")
def config_verysimple():
    filename = "tardis_configv1_verysimple.yml"
    path = os.path.abspath(os.path.join("tardis/io/tests/data/", filename))
    config = Configuration.from_yaml(path)
    return config


@pytest.fixture(scope="session")
def simulation_verysimple(config_verysimple, atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(4000)
    return sim
