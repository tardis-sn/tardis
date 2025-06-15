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


@pytest.fixture(scope="session")
def he_test_config(he_test_config_yaml):
    """Default HE workflow configuration for testing (YAML version)."""
    return he_test_config_yaml
