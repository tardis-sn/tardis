from dataclasses import dataclass

from astropy import units as u


@dataclass
class ElectronDistribution:
    """Electron temperature and density distribution.

    temperature : Quantity
        Electron temperatures.
    number_density : Quantity
        Electron number densities.
    """

    temperature: u.Quantity
    number_density: u.Quantity
