import os
import pytest
from tardis.io import config_reader
from tardis.montecarlo.base import MontecarloRunner


@pytest.fixture(scope="function")
def runner():
    config_fname = "tardis/io/tests/data/tardis_configv1_verysimply.yml"
    config = config_reader.Configuration.from_yaml(config_fname)
    runner = MontecarloRunner.from_config(config)
    return runner

