from astropy import units as u

import numpy as np

class BaseAttribute(object):
    def __init__(self):
        self._subattribute_getters={}

    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            if hasattr(self, '_value'):
                return self._value
            else:
                raise AttributeError('{0} object does not provide the '
                                     'requested attribute'.format(
                    owner.__name__))





class QuantityAttribute(BaseAttribute):
    """
    A quantity attribute is a descriptor for a model attribute that is
      represented by a quantity

    Parameter:
    ----------

    default_unit: ~str or astropy.units.Unit
        default unit that the input will be cast to
    """
    def __init__(self, default_unit):
        super(QuantityAttribute, self).__init__()
        self.default_unit = default_unit


    def __set__(self, instance, value):
        if self.default_unit is not None:
            self._value = u.Quantity(value, self.default_unit)
        else:
            self._value = value


class RadialGeometryQuantityAttribute(QuantityAttribute):

    def __init__(self, default_unit, prefix):
        super(RadialGeometryQuantityAttribute, self).__init__(default_unit)
        self.prefix = prefix
        self._subattribute_getters.update({'{0}_inner'.format(prefix):
                                    self._get_inner_value,
                               '{0}_outer'.format(prefix):
                                   self._get_outer_value})
        self._subattribute_setters = {}

    def _get_inner_value(self):
        return self._value[0:-1]

    def _get_outer_value(self):
        return self._value[1:]


class Radius1DAttribute(RadialGeometryQuantityAttribute):

    def __init__(self, default_unit=u.cm, prefix='r'):
        super(Radius1DAttribute, self).__init__(default_unit, prefix)

        self._subattribute_getters.update({'volume': self._get_volume})

        self.default_name = 'radius'

    def _get_volume(self):
        volume_0 = (4 * np.pi / 3) * self._value[0]**3
        volume = np.hstack(((volume_0.value),
                            (4 * np.pi / 3) * np.diff(self._value.value**3)))
        volume = u.Quantity(volume, volume_0.unit)

        return volume

class HomologousVelocity1D(RadialGeometryQuantityAttribute):
    def __init__(self, default_unit=u.cm / u.s, prefix='v'):
        super(HomologousVelocity1D, self).__init__(default_unit, prefix)


    def __set__(self, instance, value):
        super(HomologousVelocity1D, self).__set__(instance, value)

        try:
            time = instance.time
        except AttributeError:
            time = None

        if hasattr(instance, 'radius') and hasattr(object, 'time'):
                instance.radius = self._value * instance.time




