# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import logging
import warnings

import pyne.data

from tardis.util.colored_logger import ColoredFormatter, formatter_message

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *

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


list_of_filter = []


def logging_state(log_state, tardis_config, specific):
    """
    The configuration for the Logging status of the run() simulation object
    Called from run_tardis()
    Invoked via Function, argument 
    Proposed : Invoking via YAML & Flags

    Parameters
    ----------
    log_state: boolean or string
        False, turns the logger off for the simulation, Default
        True, allows logging for the simulation
        Log_level = ["NotSet", "Debug", "Info", "Warning", "Error", "Critical"]
        Allowed values which set the particular log level for the simulation
    """
    if specific or tardis_config["debug"]["specific_logging"]:
        specific = True
    else:
        specific = False

    if tardis_config["debug"]["logging_level"] or log_state:
        if (
            log_state.upper() == "CRITICAL"
            and tardis_config["debug"]["logging_level"]
        ):
            logging_level = tardis_config["debug"]["logging_level"]
        elif log_state:
            logging_level = log_state
            if tardis_config["debug"]["logging_level"] and log_state:
                print("Log_state & logging_level both specified")
                print("Log_state will be used for Log Level Determination\n")
        else:
            logging_level = tardis_config["debug"]["logging_level"]

    loggers = [
        logging.getLogger(name) for name in logging.root.manager.loggerDict
    ]
    if logging_level.upper() in [
        "NOTSET",
        "DEBUG",
        "INFO",
        "WARNING",
        "ERROR",
        "CRITICAL",
    ]:
        logging_levels = {
            "NOTSET": logging.NOTSET,
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
            "CRITICAL": logging.CRITICAL,
        }
        for logger in loggers:
            logger.setLevel(logging_levels[logging_level.upper()])

    if not list_of_filter == []:
        for filter in list_of_filter:
            for logger in loggers:
                logger.removeFilter(filter)

    if specific:
        filter_log = FilterLog(logging_levels[logging_level.upper()])
        list_of_filter.append(filter_log)
        for logger in loggers:
            logger.addFilter(filter_log)
    else:
        for filter in list_of_filter:
            for logger in loggers:
                logger.removeFilter(filter)


# ----------------------------------------------------------------------------
# pyne holds Python 3.7 on macOS, but refdata is pickled with protocol 5 (3.8.3)

if sys.version_info < (3, 8, 3):
    import pickle5

    sys.modules["pickle"] = pickle5
