from numba import float64
from numba.experimental import jitclass
import numpy as np
from astropy import units as u


v_x = np.linspace(-1,1, 10)*1000
v_y = np.linspace(-1,1, 10)*1000
v_z = np.linspace(-1,1, 10)*1000

r = np.dstack(np.meshgrid(v_x,v_y,v_z)).reshape(-1,3)



def cell_map(current_cell, direction):
    '''direction is 3-vector (delta_x, delta_y, delta_z)'''
    new_cell = current_cell + direction[2] + 10*direction[1] + 100 * direction[0]


class Cartesian3DGeometry:
    """
    Holds information about model geometry for Cartesian 3D models.

   Parameters
    ----------
    r_inner : astropy.units.quantity.Quantity [left back lower] [-x, -y, -z]
    r_outer : astropy.units.quantity.Quantity [right forward upper] [x, y, z]
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
        """Volume in cell"""
        return np.prod(r_outer-r_inner, axis=1)

    def to_numba(self):
        """
        Returns a new NumbaCartesian3DGeometry object

        Returns
        -------
        NumbaCartesian3DGeometry
            Numba version of Cartesian3DGeometry with properties in cgs units
        """
        return NumbaCartesian3DGeometry(
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
class NumbaCartesian3DGeometry(object):
    def __init__(self, r_inner, r_outer, v_inner, v_outer):
        """
        Cartesian 3D Geometry for the Numba mode

        Parameters

        r_inner : numpy.ndarray
        r_outer : numpy.ndarray
        v_inner : numpy.ndarray
        v_outer : numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer
