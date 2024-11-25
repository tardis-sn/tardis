from pathlib import Path

import pytest

from tardis.io.configuration.config_reader import Configuration


@pytest.fixture(scope="module")
def config(example_configuration_dir: Path):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="function")
def convergence_config(example_configuration_dir: Path):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="function")
def t_rad_strategy(convergence_config):
    return convergence_config.montecarlo.convergence_strategy.t_rad


@pytest.fixture(scope="function")
def strategy(convergence_config):
    return convergence_config.montecarlo.convergence_strategy
