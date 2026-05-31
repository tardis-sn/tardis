"""Base classes and types for TARDIS geometry."""

from enum import Enum, auto
from typing import NamedTuple

from astropy import units as u


class MeshLocation(Enum):
    """Define where properties are located on the mesh.

    Attributes
    ----------
    INTERFACE : auto
        Properties defined at cell interfaces (boundaries).
    VOLUME : auto
        Properties defined at cell volumes (centers).
    """

    INTERFACE = auto()
    VOLUME = auto()


class MeshQuantity(NamedTuple):
    """Container for a physical quantity bound to a specific location on the mesh.

    Parameters
    ----------
    data : u.Quantity
        The physical quantity data with units.
    defined_at : MeshLocation
        Location on the mesh where the quantity is defined.
    name : str, optional
        Name of the quantity, by default None.
    """

    data: u.Quantity
    defined_at: MeshLocation
    name: str = None
