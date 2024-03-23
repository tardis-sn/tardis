from numba import float64
from numba.experimental import jitclass
import numpy as np
from astropy import units as u
import warnings


class HomologousRadial1DGeometry:
    """
    Holds information about model geometry for radial 1D models.

    Parameters
    ----------
    v_inner : astropy.units.quantity.Quantity
    v_outer : astropy.units.quantity.Quantity
    v_inner_boundary : astropy.units.quantity.Quantity
    v_outer_boundary : astropy.units.quantity.Quantity

    Attributes
    ----------
    volume : astropy.units.quantity.Quantity
        Volume in each shell
    """

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

        assert "velocity" in v_inner.unit.physical_type
        assert "velocity" in v_outer.unit.physical_type

        self.v_inner = v_inner
        self.v_outer = v_outer

        # ensuring that the boundaries are within the simulation area

        if v_inner_boundary is None:
            self.v_inner_boundary = self.v_inner[0]
        elif v_inner_boundary < 0:
            warnings.warn(
                "v_inner_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self.v_inner_boundary = self.v_inner[0]
        else:
            self.v_inner_boundary = v_inner_boundary

        if v_outer_boundary is None:
            self.v_outer_boundary = self.v_outer[-1]
        elif v_outer_boundary < 0:
            warnings.warn(
                "v_outer_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self.v_outer_boundary = self.v_outer[-1]
        else:
            self.v_outer_boundary = v_outer_boundary

        assert self.v_inner_boundary < self.v_outer_boundary
        if self.v_inner_boundary < self.v_inner[0]:
            warnings.warn(
                "Requesting inner boundary below inner shell. Extrapolating the inner cell"
            )

        if self.v_outer_boundary > self.v_outer[-1]:
            warnings.warn(
                "Requesting inner boundary below inner shell. Extrapolating the inner cell"
            )

    @property
    def v_middle(self):
        return (self.v_inner + self.v_outer) / 2.0

    @property
    def v_middle_active(self):
        return self.v_middle[
            self.v_inner_boundary_index : self.v_outer_boundary_index
        ]

    @property
    def v_inner_boundary_index(self):
        return np.clip(
            np.searchsorted(self.v_inner, self.v_inner_boundary, side="right")
            - 1,
            0,
            None,
        )

    @property
    def v_outer_boundary_index(self):
        return np.clip(
            np.searchsorted(self.v_outer, self.v_outer_boundary, side="left")
            + 1,
            None,
            len(self.v_outer),
        )

    @property
    def v_inner_active(self):
        v_inner_active = self.v_inner[
            self.v_inner_boundary_index : self.v_outer_boundary_index
        ].copy()
        v_inner_active[0] = self.v_inner_boundary
        return v_inner_active

    @property
    def v_outer_active(self):
        v_outer_active = self.v_outer[
            self.v_inner_boundary_index : self.v_outer_boundary_index
        ].copy()
        v_outer_active[-1] = self.v_outer_boundary
        return v_outer_active

    @property
    def r_inner(self):
        return (self.v_inner * self.time_explosion).cgs

    @property
    def r_inner_active(self):
        return (self.v_inner_active * self.time_explosion).cgs

    @property
    def r_outer(self):
        return (self.v_outer * self.time_explosion).cgs

    @property
    def r_outer_active(self):
        return (self.v_outer_active * self.time_explosion).cgs

    @property
    def r_middle(self):
        return (self.v_middle * self.time_explosion).cgs

    @property
    def r_middle_active(self):
        return (self.v_middle_active * self.time_explosion).cgs

    @property
    def volume(self):
        """Volume in shell computed from r_outer and r_inner"""
        return (4.0 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)

    @property
    def volume_active(self):
        """Volume in shell computed from r_outer and r_inner"""
        return (
            (4.0 / 3)
            * np.pi
            * (self.r_outer_active**3 - self.r_inner_active**3)
        )

    @property
    def no_of_shells(self):
        return len(self.r_inner)

    @property
    def no_of_shells_active(self):
        return len(self.r_inner_active)

    def to_numba(self):
        """
        Returns a new NumbaRadial1DGeometry object

        Returns
        -------
        NumbaRadial1DGeometry
            Numba version of Radial1DGeometry with properties in cgs units
        """
        return NumbaRadial1DGeometry(
            self.r_inner_active.to(u.cm).value,
            self.r_outer_active.to(u.cm).value,
            self.v_inner_active.to(u.cm / u.s).value,
            self.v_outer_active.to(u.cm / u.s).value,
        )


numba_geometry_spec = [
    ("r_inner", float64[:]),
    ("r_outer", float64[:]),
    ("v_inner", float64[:]),
    ("v_outer", float64[:]),
    ("volume", float64[:]),
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
        volume : numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer
        self.volume = (4 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)
