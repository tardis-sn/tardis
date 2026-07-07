import warnings
import math

import numpy as np
from astropy import units as u
from numba import float64
from numba.experimental import jitclass

from tardis import constants as const

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


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
        return (self.v_inner_active + self.v_outer_active) / 2.0

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
            1.0 / self.time_explosion.to(u.s).value,
        )


numba_geometry_spec = [
    ("r_inner", float64[:]),
    ("r_outer", float64[:]),
    ("v_inner", float64[:]),
    ("v_outer", float64[:]),
    ("inverse_time_explosion", float64),
    ("volume", float64[:]),
]


@jitclass(numba_geometry_spec)
class NumbaRadial1DGeometry:
    def __init__(
        self, r_inner, r_outer, v_inner, v_outer, inverse_time_explosion
    ):
        """
        Radial 1D Geometry for the Numba mode

        Parameters
        ----------
        r_inner : numpy.ndarray
        r_outer : numpy.ndarray
        v_inner : numpy.ndarray
        v_outer : numpy.ndarray
        inverse_time_explosion : float
            Inverse time since explosion in 1/s.
        volume : numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer
        self.inverse_time_explosion = inverse_time_explosion
        self.volume = (4 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)

    def get_velocity(self, r: float, shell_id: int) -> float:
        """
        Calculate the homologous velocity at a given radius within a shell.

        Parameters
        ----------
        r : float
            Radius at which to calculate the velocity.
        shell_id : int
            Shell index.

        Returns
        -------
        float
            Velocity at radius ``r`` within shell ``shell_id``.
        """
        return r * self.inverse_time_explosion

    def get_doppler_factor(
        self,
        r: float,
        mu: float,
        shell_id: int,
        enable_full_relativity: bool,
    ) -> float:
        """
        Calculate the lab-to-comoving Doppler factor at radius ``r``.

        Parameters
        ----------
        r : float
            Radius at which to calculate the local homologous velocity.
        mu : float
            Packet propagation angle cosine in the lab frame.
        shell_id : int
            Shell index. Ignored for homologous geometry.
        enable_full_relativity : bool
            Flag to enable full relativistic calculations.

        Returns
        -------
        float
            Doppler factor.
        """
        beta = r * self.inverse_time_explosion / C_SPEED_OF_LIGHT
        if not enable_full_relativity:
            return 1.0 - mu * beta
        return (1.0 - mu * beta) / math.sqrt(1.0 - beta * beta)

    def get_inverse_doppler_factor(
        self,
        r: float,
        mu: float,
        shell_id: int,
        enable_full_relativity: bool,
    ) -> float:
        """
        Calculate the comoving-to-lab inverse Doppler factor at radius ``r``.

        Parameters
        ----------
        r : float
            Radius at which to calculate the local homologous velocity.
        mu : float
            Packet propagation angle cosine in the comoving frame.
        shell_id : int
            Shell index. Ignored for homologous geometry.
        enable_full_relativity : bool
            Flag to enable full relativistic calculations.

        Returns
        -------
        float
            Inverse Doppler factor.
        """
        beta = r * self.inverse_time_explosion / C_SPEED_OF_LIGHT
        if not enable_full_relativity:
            return 1.0 / (1.0 - mu * beta)
        return (1.0 + mu * beta) / math.sqrt(1.0 - beta * beta)
