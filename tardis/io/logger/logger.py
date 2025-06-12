import logging
import os
import re
import sys
from dataclasses import dataclass, field
import panel as pn
from IPython.display import display
from functools import lru_cache
import pandas as pd
from tardis.io.logger.colored_logger import ColoredFormatter

PYTHON_WARNINGS_LOGGER = logging.getLogger("py.warnings")

# During the sphinx build, we don't need the ipywidgets comms.
if 'GITHUB_ACTIONS' not in os.environ:
    pn.extension(comms="ipywidgets")
else:
    pn.extension()


def get_environment():
    """Determine the execution environment.
    
    Panel behaves differently in vscode and jupyter. In jupyter, the logger widget's display handle
    has to be updated after each log entry. In vscode, this is not needed.

    Returns
    -------
    str
        The environment name.
    """
    if any(x for x in ('VSCODE_PID', 'VSCODE') if x in os.environ):
        return 'vscode'
    return 'jupyter'

def create_logger_columns():
    """Create logger scroll columns with dynamic height.
    
    Returns
    -------
    dict
        Dictionary of scroll columns for each log level.
    """
    # Create scroll columns for each log level
    columns = {}
    
    for level in ["INFO", "WARNING/ERROR", "DEBUG"]:
        column = pn.Column(
            height=10,  # Start small
            scroll=True,
            sizing_mode='stretch_width',
            styles={
                'border': '1px solid #ddd',
                'background-color': 'white'
            }
        )
        column.max_log_entries = 1000
        column._start_height = 10
        column._max_height = 300
        columns[level] = column
    
    return columns

ENVIRONMENT = get_environment()
LOG_COLUMNS = create_logger_columns()

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

class PanelWidgetLogHandler(logging.Handler):
    """Log handler for logging to scroll columns.
    
    Parameters
    ----------
    log_columns : dict
        Dictionary of scroll columns for each log level.
    colors : dict
        Dictionary mapping log levels to display colors.
    display_widget : bool, optional
        Whether to display logs in the widget. Defaults to True.
    display_handles : dict, optional
        Dictionary of display handles for each column (jupyter environment).
    """
    def __init__(self, log_columns, colors, display_widget=True, display_handles=None):
        super().__init__()
        self.log_columns = log_columns
        self.colors = colors
        self.display_widget = display_widget
        self.display_handles = display_handles or {}
        self.environment = get_environment()
        self.stream_handler = None
        if not self.display_widget:
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setFormatter(logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"))

    def emit(self, record):
        """Process and emit a log record.
        
        Parameters
        ----------
        record : logging.LogRecord
            The log record to process and display.
        """
        if not self.display_widget:
            if self.stream_handler:
                self.stream_handler.emit(record)
            return

        if isinstance(record.msg, pd.DataFrame):
            df_str = record.msg.to_string(index=True)
            html_output = f'<pre style="white-space: pre; margin: 0;">{df_str}</pre>'
        else:
            log_entry = self.format(record)
            clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
            html_output = self._format_html_output(clean_log_entry, record)
        
        self._emit_to_columns(record.levelno, html_output)

    @staticmethod
    def _remove_ansi_escape_sequences(text):
        """Remove ANSI escape sequences from string.
        
        Parameters
        ----------
        text : str
            The text containing ANSI escape sequences.
            
        Returns
        -------
        str
            Cleaned text with ANSI escape sequences removed.
        """
        ansi_escape = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")
        return ansi_escape.sub("", text)

    def _format_html_output(self, log_entry, record):
        """Format log entry as HTML with appropriate styling.
        
        Parameters
        ----------
        log_entry : str
            The log entry text to format.
        record : logging.LogRecord
            The log record containing level information.
            
        Returns
        -------
        str
            HTML-formatted log entry.
        """
        color = self.colors.get(record.levelno, self.colors["default"])
        parts = log_entry.split(" ", 2)
        if len(parts) > 2:
            prefix, levelname, message = parts
            return f'<span>{prefix}</span> <span style="color: {color}; font-weight: bold;">{levelname}</span> {message}'
        return log_entry

    def _emit_to_columns(self, level, html_output):
        """Add log entry to appropriate scroll columns.
        
        Parameters
        ----------
        level : int
            The logging level.
        html_output : str
            The HTML-formatted log message.
        """
        level_to_output = {
            logging.WARNING: "WARNING/ERROR",
            logging.ERROR: "WARNING/ERROR",
            logging.CRITICAL: "WARNING/ERROR",
            logging.INFO: "INFO",
            logging.DEBUG: "DEBUG"
        }

        # Add to specific level column
        output_key = level_to_output.get(level)
        if output_key and output_key in self.log_columns:
            level_column = self.log_columns[output_key]
            log_pane = pn.pane.HTML(
                html_output,
                styles={
                    'font-family': 'monospace',
                    'padding': '2px 8px',
                    'margin': '0',
                    'white-space': 'pre-wrap'
                }
            )
            level_column.append(log_pane)
            
            # Dynamic height adjustment
            self._adjust_column_height(level_column)
            
            # Trim old entries
            if len(level_column) > level_column.max_log_entries:
                level_column.pop(0)
            
            # Update display handle in jupyter environment
            if (self.environment == 'jupyter' and output_key in self.display_handles 
                and self.display_handles[output_key] is not None):
                self.display_handles[output_key].update(level_column)
    
    def _adjust_column_height(self, column):
        """Dynamically adjust column height based on content.
        
        Parameters
        ----------
        column : panel.Column
            The column to adjust.
        """
        if hasattr(column, '_start_height') and hasattr(column, '_max_height'):
            # Estimate height needed: ~25px per log entry
            estimated_height = max(column._start_height, len(column) * 25)
            # Cap at max height
            new_height = min(estimated_height, column._max_height)
            column.height = new_height
    
    def close(self):
        """Close the log handler.
        """
        super().close()
        if self.stream_handler:
            self.stream_handler.flush()
            self.stream_handler.close()
            self.stream_handler = None


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

def logging_state(log_level, tardis_config, specific_log_level=None, display_logging_widget=True):
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
    tardislogger = TARDISLogger(log_columns=LOG_COLUMNS)
    tardislogger.configure_logging(log_level, tardis_config, specific_log_level)
    
    if display_logging_widget and ENVIRONMENT in ['jupyter', 'vscode']:
        if ENVIRONMENT == 'jupyter':
            # Use display handles with embed for jupyter
            display_handles = {}
            for level, column in LOG_COLUMNS.items():
                level_title = pn.pane.HTML(f"<h4 style='margin: 5px 0; color: #333;'>{level} LOGS</h4>")
                display(level_title)
                display_handles[level] = display(column.embed(), display_id=f"logger_column_{level.lower().replace('/', '_')}")
            
            # Update tardislogger with display handles
            tardislogger.display_handles = display_handles
        else:
            # Use direct display for vscode
            for level, column in LOG_COLUMNS.items():
                level_title = pn.pane.HTML(f"<h4 style='margin: 5px 0; color: #333;'>{level} LOGS</h4>")
                display(level_title)
                display(column)
    
    # Setup widget logging once after display handles are configured
    tardislogger.setup_widget_logging(display_widget=display_logging_widget)
    
    if display_logging_widget and ENVIRONMENT in ['jupyter', 'vscode']:
        return LOG_COLUMNS, tardislogger
    else:
        return None, tardislogger
