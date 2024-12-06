import logging
import os
import re
from dataclasses import dataclass, field

import panel as pn
from IPython.display import display

pn.extension()

PYTHON_WARNINGS_LOGGER = logging.getLogger("py.warnings")

def get_environment():
    """Determine the execution environment"""
    try:
        import IPython
        ipython = IPython.get_ipython()

        if ipython is None:
            return 'standard'

        # Check for VSCode specific environment variables
        if any(x for x in ('VSCODE_PID', 'VSCODE') if x in os.environ):
            return 'vscode'

        # Check if running in Jupyter notebook
        if 'IPKernelApp' in ipython.config:
            return 'jupyter'

        return 'standard'
    except:
        return 'standard'

def create_output_widget(height=300):
    return pn.Feed(
        height=height,
        styles={
            'border': '1px solid #ddd',
            'width': '100%',
            'font-family': 'monospace',
            'padding': '8px',
            'background-color': 'white'
        },
        load_buffer=1_000_000, 
        view_latest=True
    )

log_outputs = {
    "WARNING/ERROR": create_output_widget(),
    "INFO": create_output_widget(),
    "DEBUG": create_output_widget(),
    "ALL": create_output_widget(),
}

tab_order = ["ALL", "WARNING/ERROR", "INFO", "DEBUG"]
logger_widget = pn.Tabs(
    *[(title, log_outputs[title]) for title in tab_order],
    height=350,
    sizing_mode='stretch_width'
)

@dataclass
class LoggingConfig:
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
    def __init__(self, log_outputs, colors, display_widget=True):
        super().__init__()
        self.log_outputs = log_outputs
        self.colors = colors
        self.environment = get_environment()
        self.display_widget = display_widget

        # Only set up display handle for Jupyter
        if self.display_widget and self.environment == 'jupyter':
            self.display_handle = display(logger_widget, display_id=True)

    def emit(self, record):
        # Handle standard environment with simple stream output
        if not self.display_widget or self.environment == 'standard':
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"))
            stream_handler.emit(record)
            return

        # Process and emit log directly
        log_entry = self.format(record)
        clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
        html_output = self._format_html_output(clean_log_entry, record)
        self._emit_to_widget(record.levelno, html_output)

    @staticmethod
    def _remove_ansi_escape_sequences(text):
        """Remove ANSI escape sequences from string."""
        ansi_escape = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")
        return ansi_escape.sub("", text)

    def _format_html_output(self, log_entry, record):
        """Format log entry as HTML with appropriate styling."""
        color = self.colors.get(record.levelno, self.colors["default"])
        parts = log_entry.split(" ", 2)
        if len(parts) > 2:
            prefix, levelname, message = parts
            return f'<span>{prefix}</span> <span style="color: {color}; font-weight: bold;">{levelname}</span> {message}'
        return log_entry

    def _emit_to_widget(self, level, html_output):
        """Handles the widget updates using Feed component"""
        level_to_output = {
            logging.WARNING: "WARNING/ERROR",
            logging.ERROR: "WARNING/ERROR",
            logging.INFO: "INFO",
            logging.DEBUG: "DEBUG"
        }

        html_wrapped = pn.pane.HTML(f"<div style='margin: 0;'>{html_output}</div>")

        # Update specific level output
        output_key = level_to_output.get(level)
        if output_key:
            self.log_outputs[output_key].append(html_wrapped)

        # Update ALL output
        self.log_outputs["ALL"].append(html_wrapped)

        # Update Jupyter display if in jupyter environment
        if self.environment == 'jupyter':
            self.display_handle.update(logger_widget.embed())

class TARDISLogger:
    def __init__(self):
        self.config = LoggingConfig()
        self.logger = logging.getLogger("tardis")

    def configure_logging(self, log_level, tardis_config, specific_log_level=None):
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
        """
        Set up widget-based logging interface.

        Parameters
        ----------
        display_widget : bool, optional
            Whether to display the widget in GUI environments (default: True)
        """
        self.widget_handler = AsyncEmitLogHandler(
            log_outputs,
            self.config.COLORS,
            display_widget=display_widget
        )
        self.widget_handler.setFormatter(
            logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)")
        )

        self._configure_handlers()

    def _configure_handlers(self):
        """Configure logging handlers."""
        logging.captureWarnings(True)

        for logger in [self.logger, logging.getLogger()]:
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)

        self.logger.addHandler(self.widget_handler)
        PYTHON_WARNINGS_LOGGER.addHandler(self.widget_handler)

class LogFilter:
    """Filter for controlling which log levels are displayed."""

    def __init__(self, log_levels):
        self.log_levels = log_levels

    def filter(self, log_record):
        return log_record.levelno in self.log_levels

def logging_state(log_level, tardis_config, specific_log_level=None, display_logging_widget=True):
    logger = TARDISLogger()
    logger.configure_logging(log_level, tardis_config, specific_log_level)
    logger.setup_widget_logging(display_widget=display_logging_widget)

    if display_logging_widget and get_environment() == 'vscode':
        display(logger_widget)

    return logger_widget if (display_logging_widget and get_environment() in ['jupyter', 'vscode']) else None
