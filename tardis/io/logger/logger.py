import logging
import sys

from IPython.display import display

import tardis.util.panel_init as panel_init
from tardis.io.logger.colored_logger import ColoredFormatter
from tardis.io.logger.logger_widget import (
    PanelWidgetLogHandler,
    create_logger_columns,
)
from tardis.util.environment import Environment

panel_init.auto()

PYTHON_WARNINGS_LOGGER = logging.getLogger("py.warnings")

logger = logging.getLogger(__name__)

LOGGING_LEVELS = logging.getLevelNamesMapping()
LOG_COLORS = {
    logging.INFO: "#D3D3D3",
    logging.WARNING: "orange",
    logging.ERROR: "red",
    logging.CRITICAL: "orange",
    logging.DEBUG: "blue",
    "default": "black",
}
DEFAULT_LOG_LEVEL = "INFO"
DEFAULT_SPECIFIC_LOG_LEVEL = False


class TARDISLogger:
    """Main logger class for TARDIS.

    Parameters
    ----------
    log_columns : dict
        Dictionary of scroll columns for each log level.
    display_handles : dict, optional
        Dictionary of display handles for each column (jupyter environment).
    """

    def __init__(
        self,
        log_columns: dict | None = None,
        display_handles: dict | None = None,
        batch_size: int = 10,
    ) -> None:
        self.logger = logging.getLogger("tardis")
        self.log_columns = log_columns
        self.display_handles = display_handles
        self.display_ids = {}
        self.batch_size = batch_size

    def configure_logging(
        self,
        log_level: str | None,
        tardis_config: dict,
        specific_log_level: bool | None = None,
    ) -> None:
        """Configure the logging level and filtering for TARDIS loggers.

        Parameters
        ----------
        log_level : str or None
            The logging level to use (e.g., "INFO", "DEBUG"). Overrides the
            ``debug.log_level`` configuration entry when given.
        tardis_config : dict
            Configuration dictionary containing debug settings.
        specific_log_level : bool or None, optional
            Whether to enable specific log level filtering. Filtering is
            enabled if either this argument or the
            ``debug.specific_log_level`` configuration entry is True.

        Raises
        ------
        ValueError
            If an invalid log_level is provided.
        """
        debug_config = tardis_config.get("debug", {})
        logging_level = (
            log_level or debug_config.get("log_level", DEFAULT_LOG_LEVEL)
        ).upper()
        specific_log_level = specific_log_level or debug_config.get(
            "specific_log_level", DEFAULT_SPECIFIC_LOG_LEVEL
        )
        if log_level and debug_config.get("log_level"):
            self.logger.debug(
                "log_level is defined both in Functional Argument & YAML Configuration {debug section}, "
                "log_level = %s will be used for Log Level Determination",
                logging_level,
            )

        if logging_level not in LOGGING_LEVELS:
            raise ValueError(
                f"Passed Value for log_level = {logging_level} is Invalid. Must be one of the following {list(LOGGING_LEVELS)}"
            )

        level = LOGGING_LEVELS[logging_level]
        log_filter = LogFilter(level)
        tardis_loggers = [
            logging.getLogger(name)
            for name in logging.root.manager.loggerDict
            if name.startswith("tardis")
        ]
        for tardis_logger in tardis_loggers:
            tardis_logger.setLevel(level)
            for existing_filter in tardis_logger.filters[:]:
                if isinstance(existing_filter, LogFilter):
                    tardis_logger.removeFilter(existing_filter)
            if specific_log_level:
                tardis_logger.addFilter(log_filter)

    def setup_widget_logging(self, display_widget=True):
        """Set up widget-based logging interface.

        Parameters
        ----------
        display_widget : bool, optional
            Whether to display the widget in GUI environments. Default is True.
        """
        self.widget_handler = PanelWidgetLogHandler(
            log_columns=self.log_columns,
            colors=LOG_COLORS,
            display_widget=display_widget,
            display_handles=self.display_handles,
            batch_size=self.batch_size
        )
        self.widget_handler.setFormatter(
            logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)")
        )

        self._configure_handlers()

    def _configure_handlers(self):
        """Configure logging handlers.

        Removes existing handlers and adds the widget handler to the
        TARDIS logger and Python warnings logger.
        """
        logging.captureWarnings(True)

        for logger in [self.logger, logging.getLogger()]:
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)

        self.logger.addHandler(self.widget_handler)
        PYTHON_WARNINGS_LOGGER.addHandler(self.widget_handler)

    def finalize_widget_logging(self) -> None:
        """Finalize widget logging by embedding the final state.
        """
        # Embed the final state for Jupyter environments
        if (
            Environment.allows_widget_display()
            and self.display_handles
            and self.display_ids
        ):
            print("Embedding the final state for Jupyter environments")
            for level, column in self.log_columns.items():
                if (level in self.display_handles and level in self.display_ids
                    and self.display_handles[level] is not None):
                    self.display_handles[level].update(column.embed())

    def remove_widget_handler(self):
        """Remove the widget handler from the logger.
        """
        self.logger.removeHandler(self.widget_handler)
        PYTHON_WARNINGS_LOGGER.removeHandler(self.widget_handler)
        self.widget_handler.close()

    def setup_stream_handler(self):
        """Set up notebook-based logging after widget handler is removed.
        """
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(ColoredFormatter())

        self.logger.addHandler(stream_handler)
        PYTHON_WARNINGS_LOGGER.addHandler(stream_handler)

class LogFilter:
    """Filter that only passes log records at a single log level.

    Parameters
    ----------
    log_level : int
        The logging level to allow through the filter.
    """

    def __init__(self, log_level: int) -> None:
        self.log_level = log_level

    def filter(self, log_record: logging.LogRecord) -> bool:
        """Determine if a log record should be displayed.

        Parameters
        ----------
        log_record : logging.LogRecord
            The log record to evaluate.

        Returns
        -------
        bool
            True if the record's level is the filter's level, False otherwise.
        """
        return log_record.levelno == self.log_level

def logging_state(
    log_level: str | None,
    tardis_config: dict,
    specific_log_level: bool | None = None,
    display_logging_widget: bool = True,
    widget_start_height: int = 10,
    widget_max_height: int = 300,
    batch_size: int = 10,
) -> tuple[dict | None, TARDISLogger]:
    """Configure and initialize the TARDIS logging system.

    Parameters
    ----------
    log_level : str or None
        The logging level to use (e.g., "INFO", "DEBUG").
    tardis_config : dict
        Configuration dictionary containing debug settings.
    specific_log_level : bool or None, optional
        Whether to enable specific log level filtering.
    display_logging_widget : bool, optional
        Whether to display the logging widget. Default is True.
    widget_start_height : int, optional
        Starting height for widget columns. Default is 10.
    widget_max_height : int, optional
        Maximum height for widget columns. Default is 300.
    batch_size : int, optional
        Number of logs to batch before updating widget. Default is 10.

    Returns
    -------
    log_columns : dict or None
        Dictionary of log columns if the logging widget is displayed,
        otherwise None.
    tardislogger : TARDISLogger
        The configured TARDIS logger.
    """
    log_columns = create_logger_columns(start_height=widget_start_height, max_height=widget_max_height)
    tardislogger = TARDISLogger(log_columns=log_columns, batch_size=batch_size)
    tardislogger.configure_logging(log_level, tardis_config, specific_log_level)
    use_widget = display_logging_widget and Environment.allows_widget_display()

    if Environment.is_notebook() or Environment.is_sshjh() or Environment.is_sphinx():
        tardislogger.display_ids = {
            level: f"logger_column_{level.lower().replace('/', '_')}"
            for level in log_columns
        }
        tardislogger.display_handles = {
            level: display(column, display_id=tardislogger.display_ids[level])
            for level, column in log_columns.items()
        }
    elif Environment.is_vscode():
        # Use direct display for vscode (no change)
        for column in log_columns.values():
            display(column)
    elif Environment.is_terminal():
        logger.warning("Terminal environment detected, skipping logger widget")
    else:
        logger.warning("Unknown environment, skipping logger widget")

    # Setup widget logging once after display handles are configured
    tardislogger.setup_widget_logging(display_widget=display_logging_widget)

    return (log_columns if use_widget else None), tardislogger
