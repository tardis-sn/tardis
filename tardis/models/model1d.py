import numpy as np

from astropy import units as u
from tardis.models.attribute import (QuantityAttribute,
                                     RadialGeometryQuantityAttribute,
                                     Radius1DAttribute)

from tardis.models.base import BaseModel


class Radial1D(BaseModel):
    """
    Radial Model 1D
    """

    radius = Radius1DAttribute(u.cm, 'r')

    def __init__(self, radius, **kwargs):
        super(Radial1D, self).__init__()
        self._register_attribute('radius', radius)

    

class RadialHomologous1D(Radial1D):
    """
    Radial homologous 1D model.

    Parameters
    ----------

    velocity: ~astropy.units.Quantity
    
    """

    def __init__(self, velocity, t):
        self.velocity = velocity




