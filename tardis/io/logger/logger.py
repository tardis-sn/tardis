import logging
import os
import re
from dataclasses import dataclass, field
import panel as pn
from IPython.display import display
from functools import lru_cache
import pandas as pd

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

def create_output_widget(height=300):
    """Create an HTML pane for logging output.
    
    Parameters
    ----------
    height : int, optional
        The height of the pane in pixels.
        
    Returns
    -------
    panel.pane.HTML
        A Panel HTML pane configured for logging output.
    """
    return pn.pane.HTML(
        "",
        height=height,
        styles={
            'overflow-y': 'auto',
            'overflow-x': 'auto',
            'border': '1px solid #ddd',
            'width': '100%',
            'font-family': 'monospace',
            'padding': '8px',
            'background-color': 'white'
        }
    )

ENVIRONMENT = get_environment()
LOG_OUTPUTS = {
    "WARNING/ERROR": create_output_widget(),
    "INFO": create_output_widget(),
    "DEBUG": create_output_widget(),
    "ALL": create_output_widget(),
}
TAB_ORDER = ["ALL", "WARNING/ERROR", "INFO", "DEBUG"]
LOGGER_WIDGET = pn.Tabs(
    *[(title, LOG_OUTPUTS[title]) for title in TAB_ORDER],
    height=350,
    sizing_mode='stretch_width'
)

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

class AsyncEmitLogHandler(logging.Handler):
    """Asynchronous log handler for logging to a widget.
    
    Handles the formatting and display of log messages in Panel widgets.
    
    Parameters
    ----------
    log_outputs : dict
        Dictionary of Panel HTML panes for different log levels.
    colors : dict
        Dictionary mapping log levels to display colors.
    display_widget : bool, optional
        Whether to display logs in the widget. Defaults to True.
    display_handle : IPython.display.DisplayHandle, optional
        Handle for updating the display in Jupyter environment.
    logger_widget : panel.Tabs, optional
        The Panel Tabs widget containing the log outputs.
    """
    def __init__(self, log_outputs, colors, display_widget=True, display_handle=None, logger_widget=None):
        super().__init__()
        self.log_outputs = log_outputs
        self.colors = colors
        self.environment = get_environment()
        self.display_widget = display_widget
        self.display_handle = display_handle
        self.logger_widget = logger_widget

    def emit(self, record):
        """Process and emit a log record.
        
        Parameters
        ----------
        record : logging.LogRecord
            The log record to process and display.
        """
        if not self.display_widget or self.display_handle is None:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"))
            stream_handler.emit(record)
            return

        if isinstance(record.msg, pd.DataFrame):
            df_str = record.msg.to_string(index=True)
            html_output = f'<pre style="white-space: pre; margin: 0;">{df_str}</pre>'
        else:
            log_entry = self.format(record)
            clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
            html_output = self._format_html_output(clean_log_entry, record)
        self._emit_to_widget(record.levelno, html_output)

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

    def _emit_to_widget(self, level, html_output):
        """Handles the widget updates.
        
        Updates the appropriate log output widgets based on the log level
        and updates the display handle in Jupyter environments. Updates happen automatically
        in the vscode environment.
        
        Parameters
        ----------
        level : int
            The logging level (e.g., logging.INFO, logging.ERROR).
        html_output : str
            The HTML-formatted log message.
        """
        level_to_output = {
            logging.WARNING: "WARNING/ERROR",
            logging.ERROR: "WARNING/ERROR",
            logging.INFO: "INFO",
            logging.DEBUG: "DEBUG"
        }

        html_wrapped = f"<div style='margin: 0;'>{html_output}</div>"

        # Update specific level output
        output_key = level_to_output.get(level)
        if output_key:
            current = self.log_outputs[output_key].object or ""
            self.log_outputs[output_key].object = current + "\n" + html_wrapped if current else html_wrapped
        # Update ALL output
        current_all = self.log_outputs["ALL"].object or ""
        self.log_outputs["ALL"].object = current_all + "\n" + html_wrapped if current_all else html_wrapped

        # Update Jupyter display if in jupyter environment
        if self.environment == 'jupyter' and self.display_handle is not None:
            self.display_handle.update(self.logger_widget)


class TARDISLogger:
    """Main logger class for TARDIS.
    
    Handles configuration of logging levels, filters, and outputs.
    
    Parameters
    ----------
    display_handle : IPython.display.DisplayHandle, optional
        Handle for updating the display in Jupyter environment.
    logger_widget : panel.Tabs, optional
        The Panel Tabs widget containing the log outputs.
    log_outputs : dict, optional
        Dictionary of Panel HTML panes for different log levels.
    """
    def __init__(self, display_handle=None, logger_widget=None, log_outputs=None):
        self.config = LoggingConfig()
        self.logger = logging.getLogger("tardis")
        self.display_handle = display_handle
        self.logger_widget = logger_widget
        self.log_outputs = log_outputs

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
        self.widget_handler = AsyncEmitLogHandler(
            log_outputs=self.log_outputs,
            colors=self.config.COLORS,
            display_widget=display_widget,
            display_handle=self.display_handle,
            logger_widget=self.logger_widget
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
    
    Sets up the logging environment, configures log levels, and displays
    the logging widget if requested.
    
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
    panel.Tabs or None
        The logger widget if display_logging_widget is True and in a
        supported environment, otherwise None.
    """
    if display_logging_widget and ENVIRONMENT == 'jupyter':
        display_handle = display(LOGGER_WIDGET.embed(), display_id="logger_widget")
    elif display_logging_widget:
        display_handle = display(LOGGER_WIDGET, display_id="logger_widget")
    else:
        display_handle = None
    logger = TARDISLogger(display_handle=display_handle, logger_widget=LOGGER_WIDGET, log_outputs=LOG_OUTPUTS)
    logger.configure_logging(log_level, tardis_config, specific_log_level)
    logger.setup_widget_logging(display_widget=display_logging_widget)
    return LOGGER_WIDGET if (display_logging_widget and ENVIRONMENT in ['jupyter', 'vscode']) else None
 