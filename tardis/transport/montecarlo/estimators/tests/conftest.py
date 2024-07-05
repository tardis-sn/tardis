import pytest

from tardis.io.configuration.config_reader import Configuration


@pytest.fixture(scope="function")
def continuum_config(
    tardis_config_verysimple_nlte,
):
    continuum_config = Configuration.from_config_dict(
        tardis_config_verysimple_nlte
    )
    continuum_config.plasma.continuum_interaction.species = ["H I", "Ti II"]
    continuum_config.plasma.nlte_ionization_species = []

    return continuum_config
