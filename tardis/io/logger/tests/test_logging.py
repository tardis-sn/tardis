import pytest
import logging

from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.logger.logger import LOGGING_LEVELS
from tardis import run_tardis


@pytest.fixture(scope="module")
def config_verysimple(config_verysimple, atomic_data_fname):
    """
    Further simplify the `config_verysimple` fixture for testing logging.
    """
    config_verysimple.montecarlo.iterations = 2
    config_verysimple.montecarlo.no_of_packets = 400
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.num = 2000
    config_verysimple.atom_data = atomic_data_fname
    return config_verysimple


@pytest.fixture(scope="module")
def config_verysimple_logger(atomic_data_fname):
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
    )
    config["atom_data"] = atomic_data_fname
    return config


def test_logging_simulation(caplog, config_verysimple, atomic_data_fname):
    """
    Testing the logs for simulations runs
    """
    config_verysimple["atom_data"] = atomic_data_fname
    simulation = Simulation.from_config(config_verysimple)
    simulation.run()

    for record in caplog.records:
        assert record.levelno >= logging.INFO


# Testing Configuration of Logger via run_tardis() Function
@pytest.mark.parametrize(
    ["log_level", "specific_log_level"],
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
        self,
        caplog,
        config_verysimple_logger,
        log_level,
        specific_log_level,
    ):
        caplog.clear()
        run_tardis(
            config=config_verysimple_logger,
            log_level=log_level,
            specific_log_level=specific_log_level,
        )
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_config_yaml(
        self,
        caplog,
        config_verysimple_logger,
        log_level,
        specific_log_level,
    ):
        config = config_verysimple_logger
        config["debug"]["log_level"] = log_level
        config["debug"]["specific_log_level"] = specific_log_level

        caplog.clear()
        run_tardis(config=config)
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_both_specified(
        self,
        atomic_data_fname,
        caplog,
        config_verysimple_logger,
        log_level,
        specific_log_level,
    ):
        config = config_verysimple_logger
        config["debug"]["log_level"] = log_level
        config["debug"]["specific_log_level"] = specific_log_level

        caplog.clear()
        run_tardis(
            config=config,
            log_level=log_level,
            specific_log_level=specific_log_level,
        )
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]
