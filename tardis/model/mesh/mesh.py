"""Mesh classes for TARDIS geometry."""

from dataclasses import dataclass

import numpy as np
from astropy import units as u

from tardis.io.model.classic.base import MeshLocation, MeshQuantity


@dataclass
class Radial1DMesh:
    """1D radial finite-volume mesh defined by cell interface radii.

    The innermost point defines a central volume (index 0).
    Note: The central volume may be masked in downstream transport calculations.

    Parameters
    ----------
    radius : MeshQuantity
        Radii at cell interfaces with MeshLocation.INTERFACE.
    """

    radius: MeshQuantity

    @property
    def n_cells(self) -> int:
        """Total number of cells (including central volume).

        Returns
        -------
        int
            Number of cells in the mesh.
        """
        return len(self.radius.data) - 1

    @property
    def n_volumes(self) -> int:
        """Total number of volumes (same as n_cells).

        Returns
        -------
        int
            Number of volumes in the mesh.
        """
        return self.n_cells

    @property
    def cell_interface_id(self) -> np.ndarray:
        """Interface indices: [0, 1, ..., N].

        Returns
        -------
        np.ndarray
            Array of interface indices.
        """
        return np.arange(len(self.radius.data), dtype=np.int64)

    @property
    def cell_id(self) -> np.ndarray:
        """Cell/volume indices: [0, 1, ..., N-1] where 0 is the central volume.

        Returns
        -------
        np.ndarray
            Array of cell indices.
        """
        return np.arange(self.n_cells, dtype=np.int64)

    @property
    def volume_id(self) -> np.ndarray:
        """Alias for cell_id.

        Returns
        -------
        np.ndarray
            Array of volume indices.
        """
        return self.cell_id

    @classmethod
    def from_cell_interfaces(cls, radius: u.Quantity) -> "Radial1DMesh":
        """Create a Radial1DMesh from cell interface radii.

        Parameters
        ----------
        radius : u.Quantity
            Radii at cell interfaces, shape (N+1,).

        Returns
        -------
        Radial1DMesh
            The constructed mesh.
        """
        radius_mq = MeshQuantity(
            data=radius, defined_at=MeshLocation.INTERFACE, name="radius"
        )
        return cls(radius=radius_mq)


@dataclass
class HomologousRadial1DMesh:
    """IO helper for velocity-defined grids with homologous expansion.

    Used to bridge legacy/homologous formats to the spatial Radial1DMesh.

    Parameters
    ----------
    velocity : MeshQuantity
        Velocities at cell interfaces with MeshLocation.INTERFACE.
    time_explosion : u.Quantity
        Explosion time.
    """

    velocity: MeshQuantity
    time_explosion: u.Quantity

    @classmethod
    def from_velocity_interfaces(
        cls, velocity: u.Quantity, time_explosion: u.Quantity
    ) -> "HomologousRadial1DMesh":
        """Create a HomologousRadial1DMesh from velocity interfaces.

        Parameters
        ----------
        velocity : u.Quantity
            Velocities at cell interfaces, shape (N+1,).
        time_explosion : u.Quantity
            Explosion time.

        Returns
        -------
        HomologousRadial1DMesh
            The constructed homologous mesh.
        """
        velocity_mq = MeshQuantity(
            data=velocity, defined_at=MeshLocation.INTERFACE, name="velocity"
        )
        return cls(velocity=velocity_mq, time_explosion=time_explosion)

    def to_spatial_mesh(self):
        """Convert to a spatial Radial1DMesh.

        Returns
        -------
        mesh : Radial1DMesh
            The spatial mesh at t = time_explosion.
        velocity : MeshQuantity
            The velocity field.
        """
        radius = (self.velocity.data * self.time_explosion).to(u.cm)
        mesh = Radial1DMesh.from_cell_interfaces(radius)
        return mesh, self.velocity
