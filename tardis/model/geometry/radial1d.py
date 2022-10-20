from numba import float64
from numba.experimental import jitclass
import numpy as np
from astropy import units as u


class Radial1DGeometry:
    """
    Holds information about model geometry for radial 1D models.

    Parameters
    ----------
    r_inner : astropy.units.quantity.Quantity
    r_outer : astropy.units.quantity.Quantity
    v_inner : astropy.units.quantity.Quantity
    v_outer : astropy.units.quantity.Quantity

    Attributes
    ----------
    volume : astropy.units.quantity.Quantity
        Volume in each shell
    """

    def __init__(self, r_inner, r_outer, v_inner, v_outer):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer

    @property
    def volume(self):
        """Volume in shell computed from r_outer and r_inner"""
        return (4.0 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)

    def to_numba(self):
        """
        Returns a new NumbaRadial1DGeometry object

        Returns
        -------
        NumbaRadial1DGeometry
            Numba version of Radial1DGeometry with properties in cgs units
        """
        return NumbaRadial1DGeometry(
            self.r_inner.to(u.cm).value,
            self.r_outer.to(u.cm).value,
            self.v_inner.to(u.cm / u.s).value,
            self.v_outer.to(u.cm / u.s).value,
        )


numba_geometry_spec = [
    ("r_inner", float64[:]),
    ("r_outer", float64[:]),
    ("v_inner", float64[:]),
    ("v_outer", float64[:]),
]


@jitclass(numba_geometry_spec)
class NumbaRadial1DGeometry(object):
    def __init__(self, r_inner, r_outer, v_inner, v_outer):
        """
        Radial 1D Geometry for the Numba mode

        Parameters
        ----------
        r_inner : numpy.ndarray
        r_outer : numpy.ndarray
        v_inner : numpy.ndarray
        v_outer : numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer
