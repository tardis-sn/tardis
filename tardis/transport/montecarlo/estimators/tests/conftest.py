import pytest

from tardis.io.configuration.config_reader import Configuration


@pytest.fixture(scope="function")
def continuum_config(
    tardis_config_verysimple,
):
    continuum_config = Configuration.from_config_dict(
        tardis_config_verysimple
    )
    continuum_config.plasma.continuum_interaction.species = ["H I", "Ti II"]

    return continuum_config
