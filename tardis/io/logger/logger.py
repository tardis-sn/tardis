import logging
import re
from ipywidgets import Output, Tab, Layout
from IPython.display import display, HTML

LOGGING_LEVELS = {
    "NOTSET": logging.NOTSET,
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}
DEFAULT_LOG_LEVEL = "INFO"
DEFAULT_SPECIFIC_STATE = False

def logging_state(log_level, tardis_config, specific_log_level=None):
    """
    Function to set the logging configuration for the simulation output
    Called from within run_tardis()
    Configured via functional arguments passed through run_tardis() - log_level & specific_log_level
    Configured via YAML parameters under `debug` section - log_level & specific_log_level

    Parameters
    ----------
    log_level : str
        Allows input of the log level for the simulation.
        Uses Python logging framework to determine the messages that will be output.
    tardis_config : dict
        Configuration dictionary for TARDIS.
    specific_log_level : bool
        Allows setting specific logging levels. Logs of the `log_level` level would be output.
    """
    if "debug" in tardis_config:
        specific_log_level = (
            tardis_config["debug"].get("specific_log_level", specific_log_level)
        )
        logging_level = log_level or tardis_config["debug"].get("log_level", "INFO")
        if log_level and tardis_config["debug"].get("log_level"):
            print(
                "log_level is defined both in Functional Argument & YAML Configuration {debug section}"
            )
            print(
                f"log_level = {log_level.upper()} will be used for Log Level Determination\n"
            )
    else:
        tardis_config["debug"] = {}
        logging_level = log_level or DEFAULT_LOG_LEVEL
        specific_log_level = specific_log_level or DEFAULT_SPECIFIC_STATE

    logging_level = logging_level.upper()
    if logging_level not in LOGGING_LEVELS:
        raise ValueError(
            f"Passed Value for log_level = {logging_level} is Invalid. Must be one of the following {list(LOGGING_LEVELS.keys())}"
        )

    logger = logging.getLogger("tardis")
    tardis_loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict if name.startswith("tardis")]

    if logging_level in LOGGING_LEVELS:
        for logger in tardis_loggers:
            logger.setLevel(LOGGING_LEVELS[logging_level])

    if logger.filters:
        for filter in logger.filters:
            for logger in tardis_loggers:
                logger.removeFilter(filter)

    if specific_log_level:
        filter_log = FilterLog([LOGGING_LEVELS[logging_level], logging.INFO, logging.DEBUG])
        for logger in tardis_loggers:
            logger.addFilter(filter_log)
    else:
        for filter in logger.filters:
            for logger in tardis_loggers:
                logger.removeFilter(filter)

log_outputs = {
    "WARNING/ERROR": Output(layout=Layout(height='300px', overflow_y='auto')),
    "INFO": Output(layout=Layout(height='300px', overflow_y='auto')),
    "DEBUG": Output(layout=Layout(height='300px', overflow_y='auto')),
    "ALL": Output(layout=Layout(height='300px', overflow_y='auto'))
}

tab = Tab(children=[log_outputs["WARNING/ERROR"], log_outputs["INFO"], log_outputs["DEBUG"], log_outputs["ALL"]])
tab.set_title(0, "WARNING/ERROR")
tab.set_title(1, "INFO")
tab.set_title(2, "DEBUG")
tab.set_title(3, "ALL")

display(tab)

def remove_ansi_escape_sequences(text):
    """
    Remove ANSI escape sequences from a string.

    Parameters
    ----------
    text : str
        The input string containing ANSI escape sequences.

    Returns
    -------
    str
        The cleaned string without ANSI escape sequences.
    """
    ansi_escape = re.compile(r'\x1B[@-_][0-?]*[ -/]*[@-~]')
    return ansi_escape.sub('', text)

class WidgetHandler(logging.Handler):
    """
    A custom logging handler that outputs log messages to IPython widgets.

    Parameters
    ----------
    logging.Handler : class
        Inherits from the logging.Handler class.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def emit(self, record):
        """
        Emit a log record.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to be emitted.
        """
        log_entry = self.format(record)
        clean_log_entry = remove_ansi_escape_sequences(log_entry)

        if record.levelno == logging.INFO:
            color = '#D3D3D3'
        elif record.levelno == logging.WARNING:
            color = 'orange'
        elif record.levelno == logging.ERROR:
            color = 'red'
        elif record.levelno == logging.CRITICAL:
            color = 'orange'
        elif record.levelno == logging.DEBUG:
            color = 'blue'
        else:
            color = 'black'

        parts = clean_log_entry.split(' ', 2)
        if len(parts) > 2:
            prefix = parts[0]
            levelname = parts[1]
            message = parts[2]
            html_output = f'<span>{prefix}</span> <span style="color: {color}; font-weight: bold;">{levelname}</span> {message}'
        else:
            html_output = clean_log_entry

        if record.levelno in (logging.WARNING, logging.ERROR):
            with log_outputs["WARNING/ERROR"]:
                display(HTML(f"<pre style='white-space: pre-wrap; word-wrap: break-word;'>{html_output}</pre>"))
        elif record.levelno == logging.INFO:
            with log_outputs["INFO"]:
                display(HTML(f"<pre style='white-space: pre-wrap; word-wrap: break-word;'>{html_output}</pre>"))
        elif record.levelno == logging.DEBUG:
            with log_outputs["DEBUG"]:
                display(HTML(f"<pre style='white-space: pre-wrap; word-wrap: break-word;'>{html_output}</pre>"))
        with log_outputs["ALL"]:
            display(HTML(f"<pre style='white-space: pre-wrap; word-wrap: break-word;'>{html_output}</pre>"))

widget_handler = WidgetHandler()
widget_handler.setFormatter(logging.Formatter('%(name)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)'))

logging.captureWarnings(True)
logger = logging.getLogger("tardis")
logger.setLevel(logging.DEBUG)

# To fix the issue of duplicate logs
for handler in logger.handlers[:]:
    logger.removeHandler(handler)

root_logger = logging.getLogger()
for handler in root_logger.handlers[:]:
    root_logger.removeHandler(handler)

logger.addHandler(widget_handler)
logging.getLogger("py.warnings").addHandler(widget_handler)

class FilterLog(object):
    """
    Filter Log Class for Filtering Logging Output to a particular level.

    Parameters
    ----------
    log_levels : list
        List of log levels to be filtered.
    """
    def __init__(self, log_levels):
        self.log_levels = log_levels

    def filter(self, log_record):
        """
        Determine if the specified record is to be logged.

        Parameters
        ----------
        log_record : logging.LogRecord
            The log record to be filtered.

        Returns
        -------
        bool
            True if the log record's level is in the specified log_levels, False otherwise.
        """
        return log_record.levelno in self.log_levels
