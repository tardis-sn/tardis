"""Geometry classes for the TARDIS model IO.

This module provides convenience imports for all geometry-related classes.
For individual components, see:
- base.py: MeshLocation enum and MeshQuantity container
- mesh.py: Radial1DMesh and HomologousRadial1DMesh classes
"""

from tardis.io.model.classic.base import MeshLocation, MeshQuantity
from tardis.io.model.classic.mesh import HomologousRadial1DMesh, Radial1DMesh

__all__ = [
    "HomologousRadial1DMesh",
    "MeshLocation",
    "MeshQuantity",
    "Radial1DMesh",
]
