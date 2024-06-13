import glob
import os
import yaml
import pandas as pd
import pytest

# NOTE: the __githash__ variable is not available anymore
# on `version.py`. Find another way to fix it!

# from tardis import __githash__ as tardis_githash


def pytest_configure(config):
    integration_tests_configpath = config.getvalue("integration-tests")
    if integration_tests_configpath:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        config.integration_tests_config = yaml.load(
            open(integration_tests_configpath), Loader=yaml.CLoader
        )


@pytest.fixture(
    scope="class",
    params=[
        path
        for path in glob.glob(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "*")
        )
        if os.path.isdir(path)
    ],
)
def data_path(request):
    integration_tests_config = request.config.integration_tests_config
    hdf_filename = f"{os.path.basename(request.param)}.h5"
    if request.config.getoption("--generate-reference"):
        ref_path = os.path.join(
            os.path.expandvars(
                os.path.expanduser(integration_tests_config["reference"])
            ),
            tardis_githash[:7],
        )
    else:
        ref_path = os.path.join(
            os.path.expandvars(
                os.path.expanduser(integration_tests_config["reference"])
            ),
            hdf_filename,
        )

    path = {
        "config_dirpath": request.param,
        "reference_path": ref_path,
        "setup_name": hdf_filename[:-3],
        # Temporary hack for providing atom data per individual setup.
        # This url has all the atom data files hosted, for downloading.
        #        'atom_data_url': integration_tests_config['atom_data']['atom_data_url']
    }

    # For providing atom data per individual setup. Atom data can be fetched
    # from a local directory or a remote url.
    path["atom_data_path"] = os.path.expandvars(
        os.path.expanduser(integration_tests_config["atom_data_path"])
    )

    if request.config.getoption("--generate-reference") and not os.path.exists(
        path["reference_path"]
    ):
        os.makedirs(path["reference_path"])
    return path


@pytest.fixture(scope="class")
def reference(request, data_path):
    """Fixture to ingest reference data for slow test from already available
    HDF file. All data is unpacked as a collection of ``pd.Series`` and
    ``pd.DataFrames`` in a ``pd.HDFStore`` object and returned away.

    Assumed that ``data_path['reference_path']`` is a valid HDF file
    containing the reference dath for a particular setup.
    """
    # Reference data need not be loaded and provided if current test run itself
    # generates new reference data.
    if request.config.getoption("--generate-reference"):
        return
    else:
        try:
            reference = pd.HDFStore(data_path["reference_path"], "r")
        except IOError:
            raise IOError(
                f'Reference file {data_path["reference_path"]} does not exist and is needed'
                f" for the tests"
            )

        else:
            return reference
