import logging
import re
import panel as pn
from dataclasses import dataclass, field
import asyncio
import concurrent.futures
import threading
from IPython.display import display
import os

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
    LEVELS: dict[str, int] = field(default_factory=lambda: {
        "NOTSET": logging.NOTSET,
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL,
    })

    COLORS: dict[int | str, str] = field(default_factory=lambda: {
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
        print(f"[{threading.get_ident()}] Emit called at {asyncio.get_event_loop_policy().get_event_loop()}")
        
        if not self.display_widget or self.environment == 'standard':
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"))
            stream_handler.emit(record)
            return

        log_entry = self.format(record)
        clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
        html_output = self._format_html_output(clean_log_entry, record)

        print(f"[{threading.get_ident()}] Creating task with asyncio.to_thread")
        # Use asyncio.to_thread to run the widget update in a separate thread
        asyncio.create_task(
            asyncio.to_thread(
                self._emit_to_widget, 
                record.levelno, 
                html_output
            )
        )
        print(f"[{threading.get_ident()}] Task created")

    @staticmethod
    def _remove_ansi_escape_sequences(text):
        """Remove ANSI escape sequences from string."""
        ansi_escape = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")
        return ansi_escape.sub("", text)

    def _format_html_output(self, log_entry, record):
        """Format log entry as HTML with appropriate styling."""
        print(f"[{threading.get_ident()}] Formatting HTML output")
        color = self.colors.get(record.levelno, self.colors["default"])
        parts = log_entry.split(" ", 2)
        if len(parts) > 2:
            prefix, levelname, message = parts
            return f'<span>{prefix}</span> <span style="color: {color}; font-weight: bold;">{levelname}</span> {message}'
        return log_entry

    def _emit_to_widget(self, level, html_output):
        """Handles the actual widget updates"""
        print(f"[{threading.get_ident()}] _emit_to_widget started with level {level}")
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
            print(f"[{threading.get_ident()}] Updating {output_key} output")
            current = self.log_outputs[output_key].object or ""
            self.log_outputs[output_key].object = current + "\n" + html_wrapped if current else html_wrapped
            
        # Update ALL output
        print(f"[{threading.get_ident()}] Updating ALL output")
        current_all = self.log_outputs["ALL"].object or ""
        self.log_outputs["ALL"].object = current_all + "\n" + html_wrapped if current_all else html_wrapped

        if self.environment == 'jupyter':
            print(f"[{threading.get_ident()}] Updating Jupyter display")
            self.display_handle.update(logger_widget.embed())
        
        print(f"[{threading.get_ident()}] _emit_to_widget completed")


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


def logging_state(log_level, tardis_config, specific_log_level=None, display_widget=True):
    logger = TARDISLogger()
    logger.configure_logging(log_level, tardis_config, specific_log_level)
    logger.setup_widget_logging(display_widget=display_widget)
    
    if display_widget and get_environment() == 'vscode':
        display(logger_widget)
        
    return logger_widget if (display_widget and get_environment() in ['jupyter', 'vscode']) else None
