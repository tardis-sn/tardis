"""Test classes containing __slots__

Instance attributes named in ``__slots__`` can be introspected and are listed
in the Attributes section of the class documentation. Class attributes are
listed in the same section of the generated docs so docstrings should be used
to distinguish class attributes vs instance attributes. Regular instance
attributes are dynamically inserted into ``__dict__`` and cannot be reliably
introspected so they're not included in the documentation.
"""
from __future__ import annotations

__all__ = ['SlotDict', 'DerivedParam', 'DerivedSlotParam',]


class SlotDict(object):
    """
    A class that uses __slots__ and __dict__ for its attribute namespace.
    """
    __slots__ = {
        "instance_attr": "instance attribute docstring can be added here",
        "__dict__": None,   # Allows additional instance attributes to be added
    }

    class_attr = "class attribute value"
    """(class attr) this is a class attribute."""

    def __init__(self, param: str, other_param: str):
        """
        Initializes a SlotDict object.

        Parameters
        ----------
        param : str
            A parameter
        other_param : str
            Another parameter
        """

        self.instance_attr = param
        """Instance attributes declared in slots can also define their docstring
        here
        """

        if other_param is not None:
            self.other_attr = other_param
            """This instance attribute is dynamic (not declared in a slot) so
            it's not included in the docs
            """

    def my_method(self):
        """
        Prints the SlotDict parameters.
        """
        print(f"instance_attr: {self.instance_attr}")
        print(f"other_attr: {self.other_attr}")


class DerivedParam(SlotDict):
    """
    Extends SlotDict by adding an extra parameter
    """
    def __init__(self, param: str, other_param: str, extra_param: str):
        """
        Initializes a DerivedParam object.

        Parameters
        ----------
        param : str
            A parameter
        other_param : str
            Another parameter
        extra_param : str
            An extra parameter
        """
        super(DerivedParam, self).__init__(param, other_param)
        self.extra_attr = extra_param

    def derived_from_slot_class_method(self):
        """
        Prints the DerivedParam parameters.
        """
        print(f"instance_attr: {self.instance_attr}")
        print(f"other_attr: {self.other_attr}")
        print(f"extra_attr: {self.extra_attr}")


class DerivedSlotParam(SlotDict):
    """
    Extends SlotDict by adding a slot parameter
    """

    __slots__ = ('extra_attr',)

    def __init__(self, param: str, other_param: str, extra_param: str):
        """
        Initializes a DerivedSlotParam object.

        Parameters
        ----------
        param : str
            A parameter
        other_param : str
            Another parameter
        extra_param : str
            An extra parameter
        """
        super(DerivedSlotParam, self).__init__(param, other_param)
        self.extra_attr = extra_param

    def derived_from_slot_class_method(self):
        """
        Prints the DerivedSlotParam parameters.
        """
        print(f"instance_attr: {self.instance_attr}")
        print(f"other_attr: {self.other_attr}")
        print(f"extra_attr: {self.extra_attr}")
