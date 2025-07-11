import logging
import sys
from dataclasses import dataclass, field
import panel as pn
from IPython.display import display
import pandas as pd
from tardis.io.logger.colored_logger import ColoredFormatter
from tardis.util.environment import Environment
from tardis.io.logger.logger_widget import create_logger_columns, PanelWidgetLogHandler
import tardis.util.panel_init as panel_init
panel_init.auto()

PYTHON_WARNINGS_LOGGER = logging.getLogger("py.warnings")

logger = logging.getLogger(__name__)

@dataclass
class LoggingConfig:
    """Logging configuration.
    
    Attributes
    ----------
    LEVELS : dict
        The logging levels.
    COLORS : dict
        The logging colors.
    DEFAULT_LEVEL : str
        The default logging level.
    DEFAULT_SPECIFIC_STATE : bool
        The default specific log level state.
    """
    LEVELS: dict = field(default_factory=lambda: {
        "NOTSET": logging.NOTSET,
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL,
    })

    COLORS: dict = field(default_factory=lambda: {
        logging.INFO: "#D3D3D3",
        logging.WARNING: "orange",
        logging.ERROR: "red",
        logging.CRITICAL: "orange",
        logging.DEBUG: "blue",
        "default": "black",
    })

    DEFAULT_LEVEL = "INFO"
    DEFAULT_SPECIFIC_STATE = False

LOGGING_LEVELS = LoggingConfig().LEVELS

class TARDISLogger:
    """Main logger class for TARDIS.
    
    Parameters
    ----------
    log_columns : dict
        Dictionary of scroll columns for each log level.
    display_handles : dict, optional
        Dictionary of display handles for each column (jupyter environment).
    """
    def __init__(self, log_columns=None, display_handles=None):
        self.config = LoggingConfig()
        self.logger = logging.getLogger("tardis")
        self.log_columns = log_columns
        self.display_handles = display_handles
        self.display_ids = {}

    def configure_logging(self, log_level, tardis_config, specific_log_level=None):
        """Configure the logging level and filtering for TARDIS loggers.
        
        Parameters
        ----------
        log_level : str
            The logging level to use (e.g., "INFO", "DEBUG").
        tardis_config : dict
            Configuration dictionary containing debug settings.
        specific_log_level : bool, optional
            Whether to enable specific log level filtering.
            
        Raises
        ------
        ValueError
            If an invalid log_level is provided.
        """
        if "debug" in tardis_config:
            specific_log_level = tardis_config["debug"].get(
                "specific_log_level", specific_log_level
            )
            logging_level = log_level or tardis_config["debug"].get(
                "log_level", "INFO"
            )
            if log_level and tardis_config["debug"].get("log_level"):
                self.logger.debug(
                    "log_level is defined both in Functional Argument & YAML Configuration {debug section}, "
                    f"log_level = {log_level.upper()} will be used for Log Level Determination"
                )
        else:
            tardis_config["debug"] = {}
            logging_level = log_level or self.config.DEFAULT_LEVEL
            specific_log_level = specific_log_level or self.config.DEFAULT_SPECIFIC_STATE

        logging_level = logging_level.upper()
        if logging_level not in self.config.LEVELS:
            raise ValueError(
                f"Passed Value for log_level = {logging_level} is Invalid. Must be one of the following {list(self.config.LEVELS.keys())}"
            )

        logger = logging.getLogger("tardis")
        tardis_loggers = [
            logging.getLogger(name)
            for name in logging.root.manager.loggerDict
            if name.startswith("tardis")
        ]

        if logging_level in self.config.LEVELS:
            for logger in tardis_loggers:
                logger.setLevel(self.config.LEVELS[logging_level])

        if logger.filters:
            for filter in logger.filters:
                for logger in tardis_loggers:
                    logger.removeFilter(filter)

        if specific_log_level:
            filter_log = LogFilter([self.config.LEVELS[logging_level], logging.INFO, logging.DEBUG])
            for logger in tardis_loggers:
                logger.addFilter(filter_log)
        else:
            for filter in logger.filters:
                for logger in tardis_loggers:
                    logger.removeFilter(filter)

    def setup_widget_logging(self, display_widget=True):
        """Set up widget-based logging interface.

        Parameters
        ----------
        display_widget : bool, optional
            Whether to display the widget in GUI environments. Default is True.
        """
        self.widget_handler = PanelWidgetLogHandler(
            log_columns=self.log_columns,
            colors=self.config.COLORS,
            display_widget=display_widget,
            display_handles=self.display_handles
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
    
    def finalize_widget_logging(self):
        """Finalize widget logging by embedding the final state.
        """
        # Embed the final state for Jupyter environments
        if (Environment.is_notebook() and hasattr(self, 'display_handles') 
            and hasattr(self, 'display_ids') and self.display_handles and self.display_ids):
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
    """Filter for controlling which log levels are displayed.
    
    Parameters
    ----------
    log_levels : list
        List of logging levels to allow through the filter.
    """
    def __init__(self, log_levels):
        self.log_levels = log_levels

    def filter(self, log_record):
        """Determine if a log record should be displayed.
        
        Parameters
        ----------
        log_record : logging.LogRecord
            The log record to evaluate.
            
        Returns
        -------
        bool
            True if the record's level is in the allowed levels, False otherwise.
        """
        return log_record.levelno in self.log_levels

def logging_state(log_level, tardis_config, specific_log_level=None, display_logging_widget=True, widget_start_height=10, widget_max_height=300):
    """Configure and initialize the TARDIS logging system.
    
    Parameters
    ----------
    log_level : str
        The logging level to use (e.g., "INFO", "DEBUG").
    tardis_config : dict
        Configuration dictionary containing debug settings.
    specific_log_level : bool, optional
        Whether to enable specific log level filtering.
    display_logging_widget : bool, optional
        Whether to display the logging widget. Default is True.
        
    Returns
    -------
    dict
        Dictionary of log columns if display_logging_widget is True, otherwise None.
    """
    log_columns = create_logger_columns(start_height=widget_start_height, max_height=widget_max_height)
    tardislogger = TARDISLogger(log_columns=log_columns)
    tardislogger.configure_logging(log_level, tardis_config, specific_log_level)
    use_widget = display_logging_widget and (Environment.is_notebook() or Environment.is_vscode())
    
    if Environment.is_notebook():
        display_handles = {}
        display_ids = {}
        for level, column in log_columns.items():
            level_title = pn.pane.HTML(f"<h4 style='margin: 5px 0; color: #333;'>{level} LOGS</h4>")
            display(level_title)
            display_id = f"logger_column_{level.lower().replace('/', '_')}"
            display_handles[level] = display(column, display_id=display_id)
            display_ids[level] = display_id
        
        # Update tardislogger with display handles and IDs
        tardislogger.display_handles = display_handles
        tardislogger.display_ids = display_ids
    elif Environment.is_vscode():
        # Use direct display for vscode (no change)
        for level, column in log_columns.items():
            level_title = pn.pane.HTML(f"<h4 style='margin: 5px 0; color: #333;'>{level} LOGS</h4>")
            display(level_title)
            display(column)
    elif Environment.is_terminal():
        logger.warning("Terminal environment detected, skipping logger widget")
    else:
        logger.warning("Unknown environment, skipping logger widget")

    # Setup widget logging once after display handles are configured
    tardislogger.setup_widget_logging(display_widget=display_logging_widget)
    
    if use_widget:
        return log_columns, tardislogger
    else:
        return None, tardislogger
