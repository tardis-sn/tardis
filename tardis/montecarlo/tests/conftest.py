import os
import pytest
from tardis.io import config_reader
from tardis.montecarlo.base import MontecarloTransport


@pytest.fixture(scope="function")
def transport():
    config_fname = "tardis/io/tests/data/tardis_configv1_verysimply.yml"
    config = config_reader.Configuration.from_yaml(config_fname)
    transport = MontecarloTransport.from_config(config)
    return transport
