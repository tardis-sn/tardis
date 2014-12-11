from astropy import units as u

import numpy as np

class BaseAttribute(object):
    def __init__(self):
        self._subattribute_getters={}
        self._attribute_doc = None
        self.label = None

    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            if self.label in instance.__dict__:
                return instance.__dict__[self.label]
            else:
                raise AttributeError(
                    '{0} object does has no attribute {1}'.format(
                        owner.__name__, self.label))

    def __set__(self, instance, value):
        if self.label not in instance.__dict__:
            self._register_attribute(instance)

        self._check_len(instance, value)
        instance.__dict__[self.label] = value

    def _check_len(self, instance, value):
        try:
            attribute_len = len(value)
        except TypeError:
            return

        if '_attribute_len_' not in instance.__dict__:
            instance.__dict__['_attribute_len_'] = attribute_len
        else:
            if instance.__dict__['_attribute_len_'] != attribute_len:
                raise AttributeError(
                    'Setting attribute {0} with a '
                    'different length ({1}) than current '
                    'model length ({2})'.format(
                        self.label, attribute_len,
                        instance.__dict__['_attribute_len_']))



    def _register_attribute(self, instance, value):
        """
        Initialize the subattribute getter functions on the instance
        """
        instance._subattribute_getters.update(self._subattribute_getters)
        instance.model_attributes.append(self.label)


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
        self._set_quantity(instance, value)

    def _set_quantity(self, instance, value):
        if self.label not in instance.__dict__:
            self._register_attribute(instance, value)

        self._check_len(instance, value)

        if self.default_unit is not None:
            value = u.Quantity(value, self.default_unit)

        instance.__dict__[self.label] = value

class RadialGeometryQuantityAttribute(QuantityAttribute):

    def __init__(self, default_unit, prefix):
        super(RadialGeometryQuantityAttribute, self).__init__(default_unit)
        self.prefix = prefix
        self._subattribute_getters.update({'{0}_inner'.format(prefix):
                                    self._get_inner_value,
                               '{0}_outer'.format(prefix):
                                   self._get_outer_value})
        self._subattribute_setters = {}

    def _get_inner_value(self, instance):
        return instance.__dict__[self.label][0:-1]

    def _get_outer_value(self, instance):
        return instance.__dict__[self.label][1:]


class Radius1DAttribute(RadialGeometryQuantityAttribute):

    def __init__(self, default_unit=u.cm, prefix='r'):
        super(Radius1DAttribute, self).__init__(default_unit, prefix)

        self._subattribute_getters.update({'volume': self._get_volume})


    def _get_volume(self, instance):
        radius = instance.__dict__[self.label]
        volume_0 = (4 * np.pi / 3) * radius[0]**3
        volume = np.hstack(((volume_0.value),
                            (4 * np.pi / 3) * np.diff(radius.value ** 3)))
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
            pass
        else:
            instance.__class__.radius._set_quantity(instance,
                                                    self._value * time)

class HomologousTime(QuantityAttribute):

    def __init__(self, default_unit=u.s):
        super(HomologousTime, self).__init__(default_unit)


    def __set__(self, instance, value):
        super(HomologousTime, self).__set__(instance, value)

        try:
            velocity = instance.velocity
        except AttributeError:
            pass
        else:
            instance.__class__.radius._set_quantity(instance, velocity *
                                                    instance.__dict__[
                                                        self.label])

class HomologousTime0(QuantityAttribute):

    def __init__(self, default_unit=u.s):
        super(HomologousTime0, self).__init__(default_unit)

    def __set__(self, instance, value):
        super(HomologousTime0, self).__set__(instance, value)


class HomologousDensity0(QuantityAttribute):
    def __init__(self, default_unit=u.g / u.cm**3):
        super(HomologousDensity0, self).__init__(default_unit)

        self._subattribute_getters.update({'density': self._get_scaled_density})

    def _get_scaled_density(self, instance):
        return (instance.time0 / instance.time)**3 * instance.__dict__[self.label]





