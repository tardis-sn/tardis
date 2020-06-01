import logging
from functools import wraps

DEBUG_MODE = False
BUFFER = 1
ticker = 1

def log_decorator(func):
    """
    Decorator to log functions while in debug mode, i.e., when
    `debug_montecarlo` is True in the config. Works for
    `@jit'd and `@njit`'d functions, but with a significant speed
    penalty.


    TODO: in nopython mode: do I need a context manager?
    TODO: make numpy docstring.

    :param func: function to be logged.
    :return: either the function itself, if debug_mode is true, or
    """

    # to ensure that calling `help` on the decorated function still works
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if DEBUG_MODE:
            global ticker
            ticker += 1
            if ticker % BUFFER == 0:  # without a buffer, performance suffers
                logger.debug(f'Func: {func.__name__}. Input: {(args, kwargs)}')
                logger.debug(f'Output: {result}.')
        return result

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    return wrapper