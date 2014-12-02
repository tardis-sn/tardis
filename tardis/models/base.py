from astropy import units as u

from tardis.models.attribute import BaseAttribute, QuantityAttribute

class ModelMeta(type):
    """
    Metaclass for the model classes. This class ensures that the class-attribute
    `all_model_attributes` is set and one can check what kind of keyword arguments
    are allowed
    """
    def __new__(cls, name, bases, attrs):
        # find all descriptors, auto-set their labels
        attrs['all_model_attributes'] = []

        for base in bases:
            attrs['all_model_attributes'] += getattr(base,
                                                     'all_model_attributes', [])

        for n, v in attrs.items():
            if isinstance(v, BaseAttribute):
                attrs['all_model_attributes'].append(n)

        return super(ModelMeta, cls).__new__(cls, name, bases, attrs)



class BaseModel(object):
    __metaclass__ = ModelMeta

    temperature = QuantityAttribute(u.K)
    density = QuantityAttribute(u.g / u.cm**3)


