"""
A module that imports objects from the standard library.
"""
from pathlib import Path
from datetime import time


__all__ = ['Path', 'time', 'add']


def add(a, b):
    """
    Add two numbers
    """
    return a + b
