"""
A collection of useful classes and functions
"""
from collections import OrderedDict


def add(a, b):
    """
    Add two numbers
    """
    return a + b


class MixedSpam(OrderedDict):
    """
    Special spam
    """

    def eat(self, time):
        """
        Eat special spam in the required time.
        """
        pass
