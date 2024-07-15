import sys


if sys.version_info >= (3, 11):
    from importlib.resources.abc import Traversable
else:  # pragma: no cover
    from .py38 import Traversable


__all__ = ['Traversable']
