import logging
import re
from ipywidgets import Output, Tab, Layout
from IPython.display import display, HTML
from dataclasses import dataclass, field

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
        self.log_outputs = {
            "WARNING/ERROR": self._create_output_widget(),
            "INFO": self._create_output_widget(),
            "DEBUG": self._create_output_widget(),
            "ALL": self._create_output_widget(),
        }
        self.logger = logging.getLogger("tardis")
        
    def _create_output_widget(self, height="300px"):
        return Output(layout=Layout(height=height, overflow_y="auto"))
    
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
        widget_handler = LoggingHandler(self.log_outputs, self.config.COLORS)
        widget_handler.setFormatter(
            logging.Formatter("%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)")
        )
        
        self._configure_handlers(widget_handler)
        self._create_and_display_tabs()
    
    def _configure_handlers(self, widget_handler):
        """Configure logging handlers."""
        logging.captureWarnings(True)
        self.logger.setLevel(logging.DEBUG)
        
        # Clear existing handlers
        for logger in [self.logger, logging.getLogger()]:
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
        
        self.logger.addHandler(widget_handler)
        logging.getLogger("py.warnings").addHandler(widget_handler)
    
    def _create_and_display_tabs(self):
        """Create and display the logging tabs."""
        tab = Tab(children=[
            self.log_outputs[key] for key in 
            ["WARNING/ERROR", "INFO", "DEBUG", "ALL"]
        ])
        
        for i, title in enumerate(["WARNING/ERROR", "INFO", "DEBUG", "ALL"]):
            tab.set_title(i, title)
        
        display(tab)


class LoggingHandler(logging.Handler):
    def __init__(self, log_outputs, colors):
        super().__init__()
        self.log_outputs = log_outputs
        self.colors = colors
        
    def emit(self, record):
        """Emit a log record to the appropriate widget output."""
        log_entry = self.format(record)
        clean_log_entry = self._remove_ansi_escape_sequences(log_entry)
        html_output = self._format_html_output(clean_log_entry, record)
        
        self._display_log(record.levelno, html_output)
    
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
    
    def _display_log(self, level, html_output):
        """Display log message in appropriate outputs."""
        html_wrapped = f"<pre style='white-space: pre-wrap; word-wrap: break-word;'>{html_output}</pre>"
        
        level_to_output = {
            logging.WARNING: "WARNING/ERROR",
            logging.ERROR: "WARNING/ERROR",
            logging.INFO: "INFO",
            logging.DEBUG: "DEBUG"
        }
        
        output_key = level_to_output.get(level)
        if output_key:
            with self.log_outputs[output_key]:
                display(HTML(html_wrapped))
            
        # Always display in ALL output
        with self.log_outputs["ALL"]:
            display(HTML(html_wrapped))


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
