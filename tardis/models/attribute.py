from astropy import units as u

class BaseAttribute(object):
    pass


class QuantityAttribute(BaseAttribute):
    """
    A model attribute is a descriptor for
    """
    def __init__(self, default_unit=None):
        self.default_unit = default_unit

    def __get__(self, instance, owner):
        return self._value

    def __set__(self, instance, value):
        if self.default_unit is not None:
            self._value = u.Quantity(value, self.default_unit)
        else:
            self._value = value