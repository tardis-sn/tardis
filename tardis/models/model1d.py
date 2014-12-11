import numpy as np

from astropy import units as u

from tardis.models import attribute

from tardis.models.base import BaseModel


class Radial1D(BaseModel):
    """
    Radial Model 1D
    """

    radius = attribute.Radius1DAttribute(u.cm, 'r')

    def __init__(self, radius, **kwargs):
        super(Radial1D, self).__init__(**kwargs)
        self.radius = radius



    

class HomologousRadial1D(BaseModel):
    """
    Radial homologous 1D model.

    Parameters
    ----------

    velocity: ~astropy.units.Quantity
    time: ~astropy.units.Quantity
    """

    radius = attribute.Radius1DAttribute(u.cm, 'r')
    velocity = attribute.HomologousVelocity1D(u.cm / u.s, 'v')
    time = attribute.HomologousTime()
    time0 = attribute.HomologousTime0()
    density0 = attribute.HomologousDensity0()

    def __init__(self, velocity, time, time0, density0, **kwargs):
        super(HomologousRadial1D, self).__init__(**kwargs)

        self.velocity = velocity
        self.time = time
        self.time0 = time0
        self.density0 = density0

        for key, value in kwargs.items():
            setattr(self, key, value)








