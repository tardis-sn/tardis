from loguru import logger
from functools import partialmethod
import sys
import warnings


level = ""
logger.remove()
save = None
colorize = True
format = "  [ <bold><level>{level: <8}</level></bold> ][ <bold>{name}:{function}:{line}</bold> ]- {message}"

# capturing application warnings
# extra warnings come out as logs instead of normal prints
def showwarning(message, *args, **kwargs):
    name_of_class = args[0].__name__
    message = "[" + name_of_class + "]" + "[" + args[1] + "]:  " + str(message)
    logger.warning(message)


class verbosity_filter:
    """
    This is a class that helps in switching levels of the logger.
    Parameters
    ----------
    level: str
        the minimum level of the logger class
        there are 8 levels right now
        TRACE, DEBUG, INFO, SUCCESS, WARNING, TARDIS INFO, ERROR, CRITICAL
        with CRITICAL having the highest severity
        if any level is selected, then log messages from that and levels above it will be printed
        but not below it
    logger: logger object
    """

    def __init__(self, level, logger):
        self.level = level
        self.logger = logger

    def __call__(self, record):
        levelno = self.logger.level(self.level).no
        return record["level"].no >= levelno


class Rotator:
    """
    copied from: https://loguru.readthedocs.io/en/stable/resources/recipes.html#rotating-log-file-based-on-both-size-and-time
    """

    def __init__(self, *, size):
        self._size_limit = size

    def should_rotate(self, message, file):
        file.seek(0, 2)
        if file.tell() + len(message) > self._size_limit:
            return True
        return False


# Rotate file if over 100 MB
rotator = Rotator(size=1e8)


def reset_logger():
    """
    This function resets the logger by changing the level
    the level is accessible globally
    """
    filter_ = verbosity_filter(level, logger)
    logger.add(
        sys.stdout,
        filter=filter_,
        level="TRACE",
        format=format,
        colorize=colorize,
    )
    if save:
        logger.add(
            f"file_{level}.log",
            rotation=rotator.should_rotate,
            format=format,
            filter=filter_,
        )


warnings.showwarning = showwarning


#  adding custom TARDIS INFO Level, just above WARNING
# this is an additional level, which is above warning and below error
logger.level("TARDIS INFO", no=35, color="<fg #FF4500>")
logger.__class__.tardis_info = partialmethod(
    logger.__class__.log,
    "TARDIS INFO",
)


def init():
    global level
    global colorize
    global logger
    global save
    global reset_logger
    global remove_default
