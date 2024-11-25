import logging
import re
import panel as pn
from dataclasses import dataclass, field
import asyncio
import concurrent.futures
import threading
from IPython.display import display
import os
from queue import Queue
import queue
import multiprocessing

MixedFutureType = asyncio.Future | concurrent.futures.Future


# pn.extension()

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

class TardisLogger:
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


    def setup_widget_logging(self):
        """Set up widget-based logging interface."""
        # widget_handler = LoggingHandler(log_outputs, self.config.COLORS)
        widget_handler = AsyncEmitLogHandler(log_outputs, self.config.COLORS)
        widget_handler.setFormatter(
            logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)")
        )
        
        self._configure_handlers(widget_handler)
    
    def _configure_handlers(self, widget_handler):
        """Configure logging handlers."""
        logging.captureWarnings(True)
        self.logger.setLevel(logging.DEBUG)

        for logger in [self.logger, logging.getLogger()]:
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
        
        self.logger.addHandler(widget_handler)
        logging.getLogger("py.warnings").addHandler(widget_handler)
        
def is_running_in_notebook():
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            if 'IPKernelApp' in get_ipython().config:
                return True
    except ImportError:
        pass
    return False

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

class AsyncEmitLogHandler(logging.Handler):
    def __init__(self, log_outputs, colors):
        super().__init__()
        self.log_outputs = log_outputs
        self.colors = colors
        self.environment = get_environment()
        self.main_thread_id = threading.get_ident()
        self.futures = []
        
        if self.environment == 'vscode':
            self.loop = asyncio.new_event_loop()
            self.thread = threading.Thread(target=self._run_event_loop, daemon=True)
            self.thread.start()
        else:
            self.jupyter_loop = asyncio.new_event_loop()
            self.jupyter_thread = threading.Thread(target=self._run_jupyter_loop, daemon=True)
            self.jupyter_thread.start()
            self.display_handle = display(logger_widget, display_id=True)

    def _run_jupyter_loop(self):
        """Runs event loop in separate thread for Jupyter"""
        asyncio.set_event_loop(self.jupyter_loop)
        self.jupyter_loop.run_forever()
    
    def _run_event_loop(self):
        """Runs event loop in separate thread"""
        asyncio.set_event_loop(self.loop)
        self.loop.run_forever()

    def emit(self, record):
        log_entry = self.format(record)
        clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
        html_output = self._format_html_output(clean_log_entry, record)

        if self.environment == 'vscode':
            future = asyncio.run_coroutine_threadsafe(
                self._async_emit(record.levelno, html_output),
                self.loop
            )
            self.futures.append(future)
        else:
            future = asyncio.run_coroutine_threadsafe(
                self._async_emit(record.levelno, html_output),
                self.jupyter_loop
            )
            self.futures.append(future)

    def close(self):
        if self.environment == 'vscode':
            self.loop.call_soon_threadsafe(self.loop.stop)
            self.thread.join()
        else:
            self.jupyter_loop.call_soon_threadsafe(self.jupyter_loop.stop)
            self.jupyter_thread.join()
        super().close()

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

    async def _async_emit(self, level, html_output):
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

        if self.environment == 'jupyter':
            self.display_handle.update(logger_widget)


class LogFilter:
    """Filter for controlling which log levels are displayed."""
    def __init__(self, log_levels):
        self.log_levels = log_levels
        
    def filter(self, log_record):
        return log_record.levelno in self.log_levels


def logging_state(log_level, tardis_config, specific_log_level=None):
    """Configure logging state for TARDIS."""
    logger = TardisLogger()
    logger.configure_logging(log_level, tardis_config, specific_log_level)
    logger.setup_widget_logging()
    return logger_widget
