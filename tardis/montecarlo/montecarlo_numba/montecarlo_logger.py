import logging
from functools import wraps

DEBUG_MODE = True

def log_decorator(func):
    """
    Decorator to log functions while in debug mode, i.e., when
    `debug_montecarlo` is True in the config. Works for
    `@jit'd and `@njit`'d functions, but with a significant speed
    penalty.

    Questions:
        - stdout or print to file?

    TODO: How do log *args?
    TODO: How to pass kwargs to @jit, @njit?
    TODO: in nopython mode: do I need a context manager?
    TODO: Buffer?
    TODO: make numpy docstring.
    TODO: have this know debug_mode from the config.

    :param func: function to be logged.
    :return: either the function itself, if debug_mode is true, or
    """

    # to ensure that calling `help` on the decorated function still works
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug(f'Func: {func.__name__}. Input: {(args, kwargs)}')
        result = func(*args, **kwargs)
        logger.debug(f'Output: {result}.')
        return result

    if DEBUG_MODE:
        logger = logging.getLogger(__name__)

        logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            '%(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        return wrapper
    else:
        return func