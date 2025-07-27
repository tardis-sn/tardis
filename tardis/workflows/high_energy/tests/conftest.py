"""
Pytest fixtures for TARDIS HE workflow tests.
"""

import pytest
from pathlib import Path
from tardis.io.configuration.config_reader import Configuration


@pytest.fixture(scope="session")
def he_test_config_dir():
    """Get the path to the test configuration directory."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def he_test_config_yaml(he_test_config_dir):
    """Load the YAML configuration for HE workflow tests."""
    config_path = he_test_config_dir / "tardis_he_test_config.yml"
    return Configuration.from_yaml(str(config_path))


@pytest.fixture(scope="session")
def he_test_config_csvy(he_test_config_dir):
    """Load the CSVY configuration for HE workflow tests."""
    config_path = he_test_config_dir / "tardis_he_test_config_csvy.yml"
    return Configuration.from_yaml(str(config_path))


# Parameterized fixture for all config files
@pytest.fixture(
    scope="session",
    params=[
        "tardis_config_merger_2012.yml",
        "tardis_config_doubledet_2020_1a.yml",
        "tardis_config_def_n5_2014.yml",
        "tardis_config_ddt_n100.yml",
    ],
)
def he_test_configs(request, he_test_config_dir):
    """Load parameterized configurations for HE workflow tests."""
    config_path = he_test_config_dir / request.param
    config = Configuration.from_yaml(str(config_path))
    # Add config name as attribute for test identification
    config._config_name = request.param.replace(".yml", "")
    return config
