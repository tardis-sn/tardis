import logging

import pytest

from tardis.io.logger.logger import LogFilter, TARDISLogger

EMITTED_LEVELS = [
    logging.DEBUG,
    logging.INFO,
    logging.WARNING,
    logging.ERROR,
    logging.CRITICAL,
]


@pytest.fixture
def tardis_logger():
    tardis_logger = TARDISLogger()
    yield tardis_logger
    tardis_logger.configure_logging("INFO", {}, specific_log_level=False)


@pytest.fixture
def captured_levels(caplog):
    def emit_and_capture():
        caplog.clear()
        for level in EMITTED_LEVELS:
            logging.getLogger("tardis").log(level, "test message")
        return [record.levelno for record in caplog.records]

    return emit_and_capture


@pytest.mark.parametrize(
    ["log_level", "specific_log_level", "expected_levels"],
    [
        ("Info", False, EMITTED_LEVELS[1:]),
        ("INFO", True, [logging.INFO]),
        ("DEBUG", False, EMITTED_LEVELS),
        ("DEBUG", True, [logging.DEBUG]),
        ("WARNING", True, [logging.WARNING]),
        ("ERROR", False, EMITTED_LEVELS[3:]),
        ("CRITICAL", True, [logging.CRITICAL]),
    ],
)
class TestConfigureLoggingLevels:
    def test_function_arguments(
        self,
        tardis_logger,
        captured_levels,
        log_level,
        specific_log_level,
        expected_levels,
    ):
        tardis_logger.configure_logging(log_level, {}, specific_log_level)

        assert captured_levels() == expected_levels

    def test_yaml_configuration(
        self,
        tardis_logger,
        captured_levels,
        log_level,
        specific_log_level,
        expected_levels,
    ):
        tardis_config = {
            "debug": {
                "log_level": log_level,
                "specific_log_level": specific_log_level,
            }
        }

        tardis_logger.configure_logging(None, tardis_config)

        assert captured_levels() == expected_levels


def test_configure_logging_argument_overrides_yaml_log_level(
    tardis_logger, captured_levels
):
    tardis_config = {"debug": {"log_level": "DEBUG"}}

    tardis_logger.configure_logging("ERROR", tardis_config)

    assert captured_levels() == EMITTED_LEVELS[3:]


@pytest.mark.parametrize(
    ["specific_log_level_argument", "specific_log_level_config"],
    [(True, False), (False, True), (True, True)],
)
def test_configure_logging_specific_log_level_either_source(
    tardis_logger,
    captured_levels,
    specific_log_level_argument,
    specific_log_level_config,
):
    tardis_config = {
        "debug": {
            "log_level": "Warning",
            "specific_log_level": specific_log_level_config,
        }
    }

    tardis_logger.configure_logging(
        None, tardis_config, specific_log_level_argument
    )

    assert captured_levels() == [logging.WARNING]


@pytest.mark.parametrize("tardis_config", [{}, {"debug": {}}])
def test_configure_logging_defaults(
    tardis_logger, captured_levels, tardis_config
):
    tardis_logger.configure_logging(None, tardis_config)

    assert captured_levels() == EMITTED_LEVELS[1:]


def test_configure_logging_notset_defers_to_root(tardis_logger):
    tardis_logger.configure_logging("NOTSET", {})

    assert logging.getLogger("tardis").level == logging.NOTSET


def test_configure_logging_replaces_specific_filter(tardis_logger):
    tardis_logger.configure_logging("DEBUG", {}, specific_log_level=True)
    tardis_logger.configure_logging("WARNING", {}, specific_log_level=True)

    log_filters = [
        log_filter
        for log_filter in logging.getLogger("tardis").filters
        if isinstance(log_filter, LogFilter)
    ]
    assert len(log_filters) == 1


def test_configure_logging_invalid_log_level(tardis_logger):
    with pytest.raises(ValueError, match="log_level"):
        tardis_logger.configure_logging("LOUD", {})
