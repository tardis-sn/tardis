# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *  # noqa

# ----------------------------------------------------------------------------

__all__ = []

# ----------------------------------------------------------------------------

import sys
import logging
import warnings
import pyne.data
from tardis.util.colored_logger import ColoredFormatter, formatter_message

# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_config_file as yaml_load

warnings.filterwarnings("ignore", category=pyne.utils.QAWarning)

FORMAT = "[$BOLD%(name)-20s$RESET][%(levelname)-18s]  %(message)s ($BOLD%(filename)s$RESET:%(lineno)d)"
COLOR_FORMAT = formatter_message(FORMAT, True)

logging.captureWarnings(True)
logger = logging.getLogger("tardis")
logger.setLevel(logging.INFO)

console_handler = logging.StreamHandler(sys.stdout)
console_formatter = ColoredFormatter(COLOR_FORMAT)
console_handler.setFormatter(console_formatter)

logger.addHandler(console_handler)
logging.getLogger("py.warnings").addHandler(console_handler)

LOGGING_LEVELS = {
    "NOTSET": logging.NOTSET,
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}
DEFAULT_LOG_STATE = "CRITICAL"


class FilterLog(object):
    """
    Filter Log Class for Filtering Logging Output
    to a particular level

    Parameters
    ----------
    log_level : logging object
        allows to have a filter for the
        particular log_level
    """

    def __init__(self, log_level):
        self.log_level = log_level

    def filter(self, log_record):
        """
        filter() allows to set the logging level for
        all the record that are being parsed & hence remove those
        which are not of the particular level

        Parameters
        ----------
        log_record : logging.record
            which the paricular record upon which the
            filter will be applied

        Returns
        -------
        boolean : True, if the current log_record has the
            level that of the specified log_level
            False, if the current log_record doesn't have the
            same log_level as the specified one
        """
        return log_record.levelno == self.log_level


class SimulationLogger:
    """
    Implementation for the Simulation Logger
    Contains the `logging_state` method for determining the logging status for the simulation
    Contains two helper methods: add & remove log filters, nature determined by `specific`
    """

    # Setting up a list to store the Logging Filters set by logger.addFilter()
    list_of_filters = []

    @classmethod
    def add_filter_to_logger(self, filters_list, loggers_list, logging_level):
        """
        Function for adding the filters (instance of the `FilterLog` class) to the loggers present in
        the `loggers_list` list

        Parameters
        ----------
        filter_list : list
            A list containing all the filters in a particular logger
        loggers_list : list
            A list containing all the loggers that are present when we run the `run_tardis()` function
        logging_level : str
            The logging level which needs to be filtered out
        """
        filter_log = FilterLog(LOGGING_LEVELS[logging_level])
        filters_list.append(filter_log)
        for logger in loggers_list:
            logger.addFilter(filter_log)

    @classmethod
    def remove_filter_from_logger(self, filters_list, loggers_list):
        """
        Function for removing the filters from all the loggers present in the `loggers_list`

        Parameters
        ----------
        filter_list : list
            A list containing all the filters in a particular logger
        loggers_list : list
            A list containing all the loggers that are present when we run the `run_tardis()` function
        """
        for filter in filters_list:
            for logger in loggers_list:
                logger.removeFilter(filter)

    @classmethod
    def logging_state(self, tardis_config, log_state, specific):
        """
        Function to set the logging configuration for the simulation output
        Called from within run_tardis()
        Configured via functional arguments passed through run_tardis() - log_state & specific
        Configured via YAML parameters under `debug` section - logging_level & specific_logging

        Parameters
        ----------
        log_state: str
            Allows to input the log level for the simulation
            Uses Python logging framework to determine the messages that will be output
        specific: boolean
            Allows to set specific logging levels. Logs of the `log_state` level would be output.
        """

        if "debug" in tardis_config:
            specific = (
                tardis_config["debug"]["specific"]
                if specific is None
                else specific
            )

            logging_level = (
                log_state if log_state else tardis_config["debug"]["log_state"]
            )

            # Displays a message when both log_state & tardis["debug"]["log_state"] are specified
            if log_state and tardis_config["debug"]["log_state"]:
                print(
                    "log_state is defined both in Functional Argument & YAML Configuration {debug section}"
                )
                print(
                    f"log_state = {log_state.upper()} will be used for Log Level Determination\n"
                )

        else:
            if log_state:
                logging_level = log_state
            else:
                tardis_config["debug"] = {"log_state": DEFAULT_LOG_STATE}
                logging_level = tardis_config["debug"]["log_state"]

            if specific:
                specific = specific

        logging_level = logging_level.upper()
        if not logging_level in LOGGING_LEVELS:
            raise ValueError(
                f"Passed Value for log_state = {logging_level} is Invalid. Must be one of the following {list(LOGGING_LEVELS.keys())}"
            )

        loggers = [
            logging.getLogger(name) for name in logging.root.manager.loggerDict
        ]
        if logging_level in LOGGING_LEVELS:
            for logger in loggers:
                logger.setLevel(LOGGING_LEVELS[logging_level])

        if self.list_of_filters:
            self.remove_filter_from_logger(self.list_of_filters, loggers)

        if specific:
            self.add_filter_to_logger(
                self.list_of_filters,
                loggers,
                logging_level,
            )
        else:
            self.remove_filter_from_logger(self.list_of_filters, loggers)


# ----------------------------------------------------------------------------
# pyne holds Python 3.7 on macOS, but refdata is pickled with protocol 5 (3.8.3)

if sys.version_info < (3, 8, 3):
    import pickle5

    sys.modules["pickle"] = pickle5
