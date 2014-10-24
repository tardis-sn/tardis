from astropy import units as u

class BaseModel(object):

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = u.Quantity(value, self.K)

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = u.Quantity(value, u.g / u.cm**3)


