import warnings

import numpy as np
from astropy import units as u
from numba import float64
from numba.experimental import jitclass


class NonhomologousRadial1DGeometry:
    """
    Holds information about model geometry for non-homologous radial 1D models.

    Attributes
    ----------
    volume : astropy.units.quantity.Quantity
        Volume in each shell
    """

    def __init__(
        self,
        r_inner: u.Quantity,
        r_outer: u.Quantity,
        v_inner: u.Quantity,
        v_outer: u.Quantity,
        r_inner_boundary: u.Quantity,
        r_outer_boundary: u.Quantity,
        v_inner_boundary: u.Quantity,
        v_outer_boundary: u.Quantity,
    ):

        assert np.allclose(r_inner[1:], r_outer[:-1])
        assert np.allclose(v_inner[1:], v_outer[:-1])

        assert "velocity" in v_inner.unit.physical_type
        assert "velocity" in v_outer.unit.physical_type
        assert "length" in r_inner.unit.physical_type
        assert "length" in r_outer.unit.physical_type

        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer

        # ensuring that the boundaries are within the simulation area
        if r_inner_boundary is None:
            self.r_inner_boundary = self.r_inner[0]
        elif r_inner_boundary < 0:
            warnings.warn(
                "r_inner_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self.r_inner_boundary = self.r_inner[0]
        else:
            self.r_inner_boundary = r_inner_boundary

        if r_outer_boundary is None:
            self.r_outer_boundary = self.r_outer[-1]
        elif r_outer_boundary < 0:
            warnings.warn(
                "r_outer_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self.r_outer_boundary = self.r_outer[-1]
        else:
            self.r_outer_boundary = r_outer_boundary

        assert self.r_inner_boundary < self.r_outer_boundary
        if self.r_inner_boundary < self.r_inner[0]:
            warnings.warn(
                "Requesting inner boundary below inner shell. Extrapolating the inner cell"
            )

        if self.r_outer_boundary > self.r_outer[-1]:
            warnings.warn(
                "Requesting outer boundary above outer shell. Extrapolating the outer cell"
            )

        # Set velocity boundaries directly during construction so explicitly
        # supplied radius boundaries remain authoritative.
        if v_inner_boundary is None:
            self._v_inner_boundary = self.v_inner[0]
        elif v_inner_boundary < 0:
            warnings.warn(
                "v_inner_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self._v_inner_boundary = self.v_inner[0]
        else:
            self._v_inner_boundary = v_inner_boundary

        if v_outer_boundary is None:
            self._v_outer_boundary = self.v_outer[-1]
        elif v_outer_boundary < 0:
            warnings.warn(
                "v_outer_boundary < 0, assuming default value",
                DeprecationWarning,
            )
            self._v_outer_boundary = self.v_outer[-1]
        else:
            self._v_outer_boundary = v_outer_boundary

    def _radius_at_velocity(self, velocity: u.Quantity) -> u.Quantity:
        """Return the radius corresponding to a velocity boundary.

        The velocity profile is piecewise linear in radius. Boundary updates
        are therefore mapped through the first or last cell when they
        extrapolate outside the model, matching the geometry's existing
        boundary-extrapolation behavior.
        """
        velocity_boundaries = np.concatenate(
            (self.v_inner[:1], self.v_outer)
        )
        radius_boundaries = np.concatenate(
            (self.r_inner[:1], self.r_outer)
        )

        if velocity_boundaries[0] <= velocity_boundaries[-1]:
            boundary_index = np.searchsorted(
                velocity_boundaries, velocity, side="right"
            ) - 1
        else:
            boundary_index = len(velocity_boundaries) - np.searchsorted(
                velocity_boundaries[::-1], velocity, side="left"
            ) - 1

        boundary_index = np.clip(
            boundary_index, 0, len(velocity_boundaries) - 2
        )
        velocity_gradient = (
            velocity_boundaries[boundary_index + 1]
            - velocity_boundaries[boundary_index]
        ) / (
            radius_boundaries[boundary_index + 1]
            - radius_boundaries[boundary_index]
        )
        return radius_boundaries[boundary_index] + (
            velocity - velocity_boundaries[boundary_index]
        ) / velocity_gradient

    @property
    def v_inner_boundary(self) -> u.Quantity:
        """Velocity at the active inner boundary."""
        return self._v_inner_boundary

    @v_inner_boundary.setter
    def v_inner_boundary(self, velocity: u.Quantity) -> None:
        self._v_inner_boundary = velocity
        self.r_inner_boundary = self._radius_at_velocity(velocity)

    @property
    def v_outer_boundary(self) -> u.Quantity:
        """Velocity at the active outer boundary."""
        return self._v_outer_boundary

    @v_outer_boundary.setter
    def v_outer_boundary(self, velocity: u.Quantity) -> None:
        self._v_outer_boundary = velocity
        self.r_outer_boundary = self._radius_at_velocity(velocity)

    @property
    def v_middle(self):
        return (self.v_inner + self.v_outer) / 2.0

    @property
    def v_middle_active(self):
        return (self.v_inner_active + self.v_outer_active) / 2.0

    @property
    def r_middle(self):
        return (self.r_inner + self.r_outer) / 2.0

    @property
    def r_middle_active(self):
        return (self.r_inner_active + self.r_outer_active) / 2.0

    @property
    def r_inner_boundary_index(self):
        return np.clip(
            np.searchsorted(self.r_inner, self.r_inner_boundary, side="right")
            - 1,
            0,
            None,
        )

    @property
    def r_outer_boundary_index(self):
        return np.clip(
            np.searchsorted(self.r_outer, self.r_outer_boundary, side="left")
            + 1,
            None,
            len(self.r_outer),
        )

    @property
    def v_inner_boundary_index(self):
        return self.r_inner_boundary_index

    @property
    def v_outer_boundary_index(self):
        return self.r_outer_boundary_index

    @property
    def r_inner_active(self):
        r_inner_active = self.r_inner[
            self.r_inner_boundary_index : self.r_outer_boundary_index
        ].copy()
        r_inner_active[0] = self.r_inner_boundary
        return r_inner_active

    @property
    def r_outer_active(self):
        r_outer_active = self.r_outer[
            self.r_inner_boundary_index : self.r_outer_boundary_index
        ].copy()
        r_outer_active[-1] = self.r_outer_boundary
        return r_outer_active

    @property
    def v_inner_active(self):
        # Indices should be the same as those used for r boundary index
        v_inner_active = self.v_inner[
            self.r_inner_boundary_index : self.r_outer_boundary_index
        ].copy()
        v_inner_active[0] = self.v_inner_boundary
        return v_inner_active

    @property
    def v_outer_active(self):
        # Indices should be the same as those used for r boundary index
        v_outer_active = self.v_outer[
            self.r_inner_boundary_index : self.r_outer_boundary_index
        ].copy()
        v_outer_active[-1] = self.v_outer_boundary
        return v_outer_active

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
    def velocity_gradient(self):
        """Velocity gradient in each cell"""
        return (
            (self.v_outer - self.v_inner)
            / (self.r_outer - self.r_inner)
        )

    @property
    def velocity_gradient_active(self):
        """Velocity gradient in each cell"""
        return (
            (self.v_outer_active - self.v_inner_active)
            / (self.r_outer_active - self.r_inner_active)
        )

    def get_velocity(self, r: u.Quantity, shell_id: int):
        """
        Calculate the velocity at a given radius within a shell, assuming
        a piece-wise linear velocity with radius.

        Parameters
        ----------
        r : astropy.units.quantity.Quantity
            Radius at which to calculate the velocity
        shell_id : int
            Shell index

        Returns
        -------
        astropy.units.quantity.Quantity
            Velocity at radius r within shell shell_id

        Examples
        --------
            geometry.get_velocity(rpacket.r, rpacket.current_shell_id)
        """
        return self.v_inner[shell_id] + self.velocity_gradient[shell_id] * (r - self.r_inner[shell_id])

    @property
    def no_of_shells(self):
        return len(self.r_inner)

    @property
    def no_of_shells_active(self):
        return len(self.r_inner_active)

    def to_numba(self):
        """
        Returns a new NumbaNonhomologousRadial1DGeometry object

        Returns
        -------
        NumbaNonhomologousRadial1DGeometry
            Numba version of Radial1DGeometry with properties in cgs units
        """
        return NumbaNonhomologousRadial1DGeometry(
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
    ("velocity_gradient", float64[:]),
    ("volume", float64[:]),
]


@jitclass(numba_geometry_spec)
class NumbaNonhomologousRadial1DGeometry:
    def __init__(self, r_inner, r_outer, v_inner, v_outer):
        """
        Radial 1D Geometry for the Numba mode

        Parameters
        ----------
        r_inner : numpy.ndarray
        r_outer : numpy.ndarray
        v_inner : numpy.ndarray
        v_outer : numpy.ndarray
        velocity_gradient : numpy.ndarray
        volume : numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.v_outer = v_outer
        self.velocity_gradient = (self.v_outer - self.v_inner) / (self.r_outer - self.r_inner)
        self.volume = (4 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)

    def get_velocity(self, r: float, shell_id: int):
        """
        Calculate the velocity at a given radius within a shell, assuming
        a piece-wise linear velocity with radius.

        Parameters
        ----------
        r : float
            Radius at which to calculate the velocity
        shell_id : int
            Shell index

        Returns
        -------
        float
            Velocity at radius r within shell shell_id
        """
        return self.v_inner[shell_id] + self.velocity_gradient[shell_id] * (r - self.r_inner[shell_id])
