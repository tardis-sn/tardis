DEBUG_MODE = False

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
    if DEBUG_MODE:
        logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            '%(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        def wrapper(*args, **kwargs):
            logger.debug(f'Func: {func}. Input: {(args, kwargs)}')
            result = func(*args)
            logger.debug(f'Output: {result}.')
            return result

        return wrapper
    else:
        return func