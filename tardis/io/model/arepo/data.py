from dataclasses import dataclass, field

import numpy as np
from astropy import units as u


@dataclass
class ArepoData:
    """
    Data structure for Arepo snapshot data.

    Parameters
    ----------
    time : astropy.units.Quantity
        Time of the snapshot.
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

    time: u.Quantity
    position: np.ndarray
    velocities: np.ndarray
    densities: np.ndarray
    mass: np.ndarray
    isotope_dict: dict = field(default_factory=dict)

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
