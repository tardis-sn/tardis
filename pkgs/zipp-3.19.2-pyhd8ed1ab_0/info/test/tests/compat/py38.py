import sys


if sys.version_info >= (3, 9):
    from importlib.abc import Traversable
else:  # pragma: no cover
    from importlib_resources.abc import Traversable


__all__ = ['Traversable']
