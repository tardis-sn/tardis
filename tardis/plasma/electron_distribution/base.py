from dataclasses import dataclass

from astropy import units as u


@dataclass
class ElectronDistribution:
    """Electron temperature and density distribution.

    temperature : Quantity
        Electron temperatures in Kelvin.
    number_density : Quantity
        Electron number densities in g/cm^3.
    """

    temperature: u.Quantity
    number_density: u.Quantity
