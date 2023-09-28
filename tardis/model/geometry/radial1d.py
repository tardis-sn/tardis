from numba import float64
from numba.experimental import jitclass
import numpy as np
from astropy import units as u


class HomologousRadial1DGeometry:
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

    DEFAULT_VELOCITY_UNIT = u.km / u.s
    DEFAULT_DISTANCE_UNIT = u.km

    def __init__(
        self,
        v_inner,
        v_outer,
        v_inner_boundary,
        v_outer_boundary,
        time_explosion,
    ):
        self.time_explosion = time_explosion

        # ensuring that the cells are continuous
        assert np.allclose(v_inner[1:], v_outer[:-1])

        self.v_inner = v_inner.to(self.DEFAULT_VELOCITY_UNIT)
        self.v_outer = v_outer.to(self.DEFAULT_VELOCITY_UNIT)

        # ensuring that the boundaries are within the simulation area
        assert v_inner_boundary < v_outer_boundary
        assert v_inner_boundary >= self.v_inner[0] #TBD - we could just extrapolate
        assert v_outer_boundary <= self.v_outer[-1] 

        self.v_inner_boundary = v_inner_boundary
        self.v_outer_boundary = v_outer_boundary

    @property
    def v_inner_boundary_index(self):
        # TODO potentially rename to v_inner_active_index??

        # fix to ensure that we get the same index if we are on the shell
        # boundary vs slightly above

        return np.max(
            [
                np.searchsorted(
                    self.v_inner, self.v_inner_boundary, side="left"
                )
                - 1,
                0
            ]
        )

    @property
    def v_outer_boundary_index(self):
        return np.searchsorted(self.v_outer, self.v_outer_boundary, side='right')

    @property
    def r_inner(self):
        return (self.v_inner * self.time_explosion).to(
            self.DEFAULT_DISTANCE_UNIT
        )

    @property
    def r_outer(self):
        return (self.v_inner * self.time_explosion).to(
            self.DEFAULT_DISTANCE_UNIT
        )

    @property
    def volume(self):
        """Volume in shell computed from r_outer and r_inner"""
        return (4.0 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)

    @property
    def no_of_shells(self):
        return len(self.r_inner)

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
