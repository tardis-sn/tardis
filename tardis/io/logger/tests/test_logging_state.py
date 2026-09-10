from unittest.mock import Mock

import pytest

from tardis.io.logger import logger


def test_logging_state_uses_stream_handler_without_widget_display(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tardis_logger = Mock()
    logger_factory = Mock(return_value=tardis_logger)
    create_logger_columns = Mock()
    monkeypatch.setattr(logger, "TARDISLogger", logger_factory)
    monkeypatch.setattr(logger, "create_logger_columns", create_logger_columns)
    monkeypatch.setattr(logger.Environment, "allows_widget_display", lambda: False)
    monkeypatch.setattr(logger.Environment, "is_notebook", lambda: False)
    monkeypatch.setattr(logger.Environment, "is_sshjh", lambda: False)
    monkeypatch.setattr(logger.Environment, "is_sphinx", lambda: False)
    monkeypatch.setattr(logger.Environment, "is_vscode", lambda: False)
    monkeypatch.setattr(logger.Environment, "is_terminal", lambda: True)

    log_columns, actual_logger = logger.logging_state(
        "INFO", {}, display_logging_widget=True
    )

    assert log_columns is None
    assert actual_logger is tardis_logger
    create_logger_columns.assert_not_called()
    tardis_logger.setup_stream_handler.assert_called_once_with()
    tardis_logger.setup_widget_logging.assert_not_called()


def test_logging_state_does_not_create_widget_when_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tardis_logger = Mock()
    logger_factory = Mock(return_value=tardis_logger)
    create_logger_columns = Mock()
    monkeypatch.setattr(logger, "TARDISLogger", logger_factory)
    monkeypatch.setattr(logger, "create_logger_columns", create_logger_columns)
    monkeypatch.setattr(logger.Environment, "allows_widget_display", lambda: True)

    log_columns, actual_logger = logger.logging_state(
        "INFO", {}, display_logging_widget=False
    )

    assert log_columns is None
    assert actual_logger is tardis_logger
    create_logger_columns.assert_not_called()
    tardis_logger.setup_stream_handler.assert_called_once_with()
    tardis_logger.setup_widget_logging.assert_not_called()
