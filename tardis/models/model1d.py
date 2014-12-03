import numpy as np

from astropy import units as u
from tardis.models.attribute import (QuantityAttribute,
                                     RadialGeometryQuantityAttribute,
                                     Radius1DAttribute)

from tardis.models.base import BaseModel


class RadialModel1D(BaseModel):
    """
    Radial Model 1D
    """

    radius = RadialGeometryQuantityAttribute(u.cm, 'r')

    def __init__(self, radius, **kwargs):
        super(RadialModel1D, self).__init__()
        self._register_attribute('radius', radius)


    @property
    def mass(self):
        return self.volume * self.density


class RadialHomologousModel1D(RadialModel1D):
    """
    Radial homologous 1D model.

    Parameters
    ----------

    velocity: ~astropy.units.Quantity
    
    """

    def __init__(self, velocity, temperature, density, abundances, t_0, t_model):
        self.velocity = velocity
        self.temperature = temperature
        self.density

    @property
    def radius(self):
        return self._radius

    @property
    def velocity(self):
        return self._velocity

    @property
    def t_0(self):
        return self.t_0

    @radius.setter
    def radius(self, value):
        self._radius = u.Quantity(value, unit=u.cm)
        self._velocity = self._radius / self._t_0

    @velocity.setter
    def velocity(self, value):
        self._velocity = u.Quantity(value, unit=u.cm / u.s)
        self._radius = self.t_0 * self._velocity

    @t_0.setter
    def t_0(self, value):
        self._t_0 = u.Quantity(value, u.s)
        self._radius = self.velocity * self._t_0



