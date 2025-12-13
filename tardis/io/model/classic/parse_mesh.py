"""Parse classic TARDIS configuration to mesh objects.

This module provides parsers for the "specific" structure type only.
File-based configurations should continue using the existing geometry parsers.
"""

from tardis.model.mesh import HomologousRadial1DMesh
from tardis.util.base import quantity_linspace


def parse_homologous_mesh_from_config(config, time_explosion):
    """Parse a HomologousRadial1DMesh from a TARDIS config with 'specific' structure.

    This function only works for structure.type == "specific". For file-based
    configurations, use the existing parse_geometry_from_config instead.

    Parameters
    ----------
    config : object
        Configuration object with structure.type == "specific".
    time_explosion : u.Quantity
        The time of the explosion.

    Returns
    -------
    HomologousRadial1DMesh
        The parsed homologous mesh with velocity interfaces.

    Raises
    ------
    ValueError
        If structure.type is not "specific".
    """
    structure_config = config.model.structure

    if structure_config.type != "specific":
        raise ValueError(
            f"parse_homologous_mesh_from_config only supports structure.type='specific', "
            f"got '{structure_config.type}'. Use parse_geometry_from_config for file-based configs."
        )

    velocity = quantity_linspace(
        structure_config.velocity.start,
        structure_config.velocity.stop,
        structure_config.velocity.num + 1,
    ).cgs

    return HomologousRadial1DMesh.from_velocity_interfaces(
        velocity, time_explosion
    )
