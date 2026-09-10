from pathlib import Path
from unittest.mock import patch

from tardis.io.logger.logger import (
    PYTHON_WARNINGS_LOGGER,
    TARDISLogger,
    logging_state,
)


def test_file_handler_writes_plain_utf8_logs(tmp_path: Path) -> None:
    """Write TARDIS log records to an uncolored UTF-8 file."""
    tardis_logger = TARDISLogger()
    log_path = tmp_path / "tardis.log"

    file_handler = tardis_logger.setup_file_handler(log_path)
    tardis_logger.logger.warning("log message written to file")
    file_handler.close()

    tardis_logger.logger.removeHandler(file_handler)
    PYTHON_WARNINGS_LOGGER.removeHandler(file_handler)

    log_contents = log_path.read_text(encoding="utf-8")
    assert "log message written to file" in log_contents
    assert "\x1b" not in log_contents


def test_logging_configuration_writes_to_log_file(tmp_path: Path) -> None:
    """Configure file logging through the debug configuration section."""
    log_path = tmp_path / "configured.log"
    configuration = {"debug": {"log_file": str(log_path)}}

    with (
        patch.object(TARDISLogger, "setup_widget_logging"),
        patch(
            "tardis.io.logger.logger.Environment.allows_widget_display",
            return_value=False,
        ),
        patch(
            "tardis.io.logger.logger.Environment.is_notebook",
            return_value=False,
        ),
        patch(
            "tardis.io.logger.logger.Environment.is_sshjh", return_value=False
        ),
        patch(
            "tardis.io.logger.logger.Environment.is_sphinx", return_value=False
        ),
        patch(
            "tardis.io.logger.logger.Environment.is_vscode", return_value=False
        ),
        patch(
            "tardis.io.logger.logger.Environment.is_terminal",
            return_value=False,
        ),
    ):
        _, tardis_logger = logging_state(None, configuration)

    file_handler = next(
        handler
        for handler in tardis_logger.logger.handlers
        if getattr(handler, "baseFilename", None) == str(log_path)
    )
    tardis_logger.logger.warning("configured file logging")
    file_handler.close()

    tardis_logger.logger.removeHandler(file_handler)
    PYTHON_WARNINGS_LOGGER.removeHandler(file_handler)

    assert "configured file logging" in log_path.read_text(encoding="utf-8")
