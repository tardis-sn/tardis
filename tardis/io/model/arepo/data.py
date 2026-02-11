from dataclasses import dataclass, field

import numpy as np
from astropy import units as u


@dataclass
class ArepoData:
    """
    Data structure for Arepo snapshot data.

    Parameters
    ----------
    time : astropy.units.Quantity or float
        Time of the snapshot. If a float is provided, it is assumed to be in seconds.
    pos : numpy.ndarray
        Position array in Cartesian coordinates (3, N) in cm.
    vel : numpy.ndarray
        Velocity array in Cartesian coordinates (3, N) in cm/s.
    rho : numpy.ndarray
        Density array in g/cm^3.
    mass : numpy.ndarray
        Mass array in g.
    nuc_dict : dict
        Dictionary of nuclear mass fractions keyed by species name.
    """

    time: u.Quantity | float
    position: np.ndarray
    velocities: np.ndarray
    densities: np.ndarray
    mass: np.ndarray
    isotope_dict: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Ensure time has astropy units."""
        self.time = u.Quantity(self.time, u.s)

    @property
    def volume(self) -> np.ndarray:
        """
        Calculate volume from mass and density.

        Returns
        -------
        numpy.ndarray
            Volume array in cm^3.
        """
        return self.mass / self.densities

    @property
    def species(self) -> list[str]:
        """
        Get list of species names.

        Returns
        -------
        list of str
            Species names from nuc_dict keys.
        """
        return list(self.isotope_dict.keys())
