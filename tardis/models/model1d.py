import numpy as np

from astropy import units as u
from tardis.models.attribute import (QuantityAttribute,
                                     RadialGeometryQuantityAttribute,
                                     Radius1DAttribute, HomologousVelocity1D,
                                     HomologousTime)

from tardis.models.base import BaseModel


class Radial1D(BaseModel):
    """
    Radial Model 1D
    """

    radius = Radius1DAttribute(u.cm, 'r')

    def __init__(self, radius, **kwargs):
        super(Radial1D, self).__init__(**kwargs)
        self._register_attribute('radius', radius)



    

class RadialHomologous1D(BaseModel):
    """
    Radial homologous 1D model.

    Parameters
    ----------

    velocity: ~astropy.units.Quantity
    time: ~astropy.units.Quantity
    """

    radius = Radius1DAttribute(u.cm, 'r')
    velocity = HomologousVelocity1D(u.cm / u.s, 'v')
    time = HomologousTime(u.s)

    def __init__(self, velocity, time, **kwargs):
        super(RadialHomologous1D, self).__init__(**kwargs)
        self._register_attribute('velocity', velocity)
        self._register_attribute('time', time)
        self._register_attribute('radius', self.radius)





