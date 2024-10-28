from astropy import units as u


class ElectronDistribution:
    def __init__(
        self,
        temperature: u.Quantity,
        number_density: u.Quantity,
    ):
        """Initializes the electron distribution.

        Parameters
        ----------
        temperatures : float
            Electron temperature.
        number_densities : float
            Electron number density.
        """
        self.temperature = temperature
        self.number_density = number_density
