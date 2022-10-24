"""Configure Test Suite.

This file is used to configure the behavior of pytest when using the Astropy
test infrastructure. It needs to live inside the package in order for it to
get picked up when running the tests inside an interpreter using
packagename.test

"""

import os

from astropy.version import version as astropy_version

# For Astropy 3.0 and later, we can use the standalone pytest plugin
if astropy_version < "3.0":
    from astropy.tests.pytest_plugins import *  # noqa

    del pytest_report_header
    ASTROPY_HEADER = True
else:
    try:
        from pytest_astropy_header.display import (
            PYTEST_HEADER_MODULES,
            TESTED_VERSIONS,
        )

        ASTROPY_HEADER = True
    except ImportError:
        ASTROPY_HEADER = False


def pytest_configure(config):
    """Configure Pytest with Astropy.

    Parameters
    ----------
    config : pytest configuration

    """
    if ASTROPY_HEADER:

        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        PYTEST_HEADER_MODULES.pop("Pandas", None)
        PYTEST_HEADER_MODULES["scikit-image"] = "skimage"

        from . import __version__

        packagename = os.path.basename(os.path.dirname(__file__))
        TESTED_VERSIONS[packagename] = __version__


# Uncomment the last two lines in this block to treat all DeprecationWarnings as
# exceptions. For Astropy v2.0 or later, there are 2 additional keywords,
# as follow (although default should work for most cases).
# To ignore some packages that produce deprecation warnings on import
# (in addition to 'compiler', 'scipy', 'pygments', 'ipykernel', and
# 'setuptools'), add:
#     modules_to_ignore_on_import=['module_1', 'module_2']
# To ignore some specific deprecation warning messages for Python version
# MAJOR.MINOR or later, add:
#     warnings_to_ignore_by_pyver={(MAJOR, MINOR): ['Message to ignore']}
# from astropy.tests.helper import enable_deprecations_as_exceptions  # noqa
# enable_deprecations_as_exceptions()

# -------------------------------------------------------------------------
# Here the TARDIS testing stuff begins
# -------------------------------------------------------------------------

import pytest
import pandas as pd
from tardis.io.util import yaml_load_file, YAMLLoader
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation


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


@pytest.fixture(scope="function")
def tardis_config_verysimple():
    return yaml_load_file(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml", YAMLLoader
    )


@pytest.fixture(scope="function")
def tardis_config_verysimple_nlte():
    return yaml_load_file(
        "tardis/io/tests/data/tardis_configv1_nlte.yml", YAMLLoader
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


@pytest.fixture(scope="function")
def config_montecarlo_1e5_verysimple():
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
