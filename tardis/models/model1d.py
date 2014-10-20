import numpy as np
from astropy import units as u

from tardis.models.base import BaseModel


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

    @property
    def volume(self):
        volume_0 = (4 * np.pi / 3) * self.radius[0]**3
        volume = np.hstack(((volume_0.value),
                            (4 * np.pi / 3) * np.diff(self.radius.value**3)))
        volume = u.Quantity(volume, volume_0.unit)

        return volume

    @property
    def mass(self):
        return self.volume * self.density


class RadialHomologousModel1D(RadialModel1D):

    def __init__(self, velocity, temperature, density, abundances, t_0, n_electron=None):
        pass

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



