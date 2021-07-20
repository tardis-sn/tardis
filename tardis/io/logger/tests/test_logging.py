import pytest
import logging

from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.logger.logger import LOGGING_LEVELS
from tardis import run_tardis


def test_logging_simulation(atomic_data_fname, caplog):
    """
    Testing the logs for simulations runs
    """
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )
    config["atom_data"] = atomic_data_fname

    simulation = Simulation.from_config(config)

    simulation.run()

    for record in caplog.records:
        assert record.levelno >= logging.INFO


# Testing Configuration of Logger via run_tardis() Function
@pytest.mark.parametrize(
    ["log_level", "specific"],
    [
        ("Info", False),
        ("INFO", False),
        ("INFO", True),
        ("DEBUG", False),
        ("DEBUG", True),
        ("WARNING", True),
        ("ERROR", True),
        ("CRITICAL", True),
        ("NOTSET", False),
    ],
)
class TestSimulationLogging:
    """
    Class implemented for testing the logging configuration available via run_tardis()
    Tests Functional Arguments : log_level & specific
    Tests YAML Parameters : logging_level & specific_logging
    """

    def test_logging_config(
        self, atomic_data_fname, caplog, log_level, specific
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname

        caplog.clear()
        run_tardis(config=config, log_level=log_level, specific=specific)
        for record in caplog.records:
            if specific == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_config_yaml(
        self, atomic_data_fname, caplog, log_level, specific
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname
        config["debug"]["log_level"] = log_level
        config["debug"]["specific"] = specific

        caplog.clear()
        run_tardis(config=config)
        for record in caplog.records:
            if specific == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_both_specified(
        self, atomic_data_fname, caplog, log_level, specific
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname
        config["debug"]["log_level"] = log_level
        config["debug"]["specific"] = specific

        caplog.clear()
        run_tardis(config=config, log_level=log_level, specific=specific)
        for record in caplog.records:
            if specific == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]
