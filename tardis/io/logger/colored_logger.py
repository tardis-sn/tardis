import logging

"""
Code for Custom Logger Classes (ColoredFormatter and ColorLogger) and its helper function
(formatter_message) is used from this thread
http://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
"""


def formatter_message(message, use_color=True):
    """
    Helper Function used for Coloring Log Output
    """
    # These are the sequences need to get colored ouput
    RESET_SEQ = "\033[0m"
    BOLD_SEQ = "\033[1m"
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace(
            "$BOLD", BOLD_SEQ
        )
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message


class ColoredFormatter(logging.Formatter):
    """
    Custom logger class for changing levels color
    """

    non_debug = "[$BOLD{name:20s}$RESET][{levelname:18s}]  \n\t{message:s} ($BOLD{filename:s}$RESET:{lineno:d})"
    debug = "[$BOLD{name:20s}$RESET][{levelname:18s}]  {message:s} ($BOLD{filename:s}$RESET:{lineno:d})"

    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg, style="{")
        self.use_color = use_color

    def format(self, record):
        non_debug = formatter_message(ColoredFormatter.non_debug, True)
        debug = formatter_message(ColoredFormatter.debug, True)
        COLOR_SEQ = "\033[1;%dm"
        RESET_SEQ = "\033[0m"
        BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
        COLORS = {
            "WARNING": YELLOW,
            "INFO": WHITE,
            "DEBUG": BLUE,
            "CRITICAL": YELLOW,
            "ERROR": RED,
        }
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = (
                COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            )
            record.levelname = levelname_color
        if record.levelno == logging.DEBUG:
            self._style._fmt = debug
        else:
            self._style._fmt = non_debug
        return logging.Formatter.format(self, record)


class ColoredLogger(logging.Logger):
    """
    Custom logger class with multiple destinations
    """

    FORMAT = "[$BOLD{name:20s}$RESET][{levelname:18s}]  \n\t{message:s} ($BOLD{filename:s}$RESET:{lineno:d})"
    COLOR_FORMAT = formatter_message(FORMAT, True)

    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return
