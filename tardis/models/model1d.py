from tardis.models.base import BaseModel
from astropy import units as u

class RadialModel1D(BaseModel):
    """
    Radial Model 1D

    """

    def __init__(self, radius, density, abundances, n_electron=None):
        self.radius = radius
        self.density = density

    @property
    def radius(self):
        return self._radius


    @radius.setter
    def radius(self, value):
        self._radius = u.Quantity(value, u.cm)

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = u.Quantity(value, u.g / u.cm**3)


