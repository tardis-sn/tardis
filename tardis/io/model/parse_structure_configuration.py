import logging
import os

from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_config_v1_density,
)
from tardis.io.model.readers.base import read_density_file
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.util.base import quantity_linspace

logger = logging.getLogger(__name__)

def parse_structure_config(config, time_explosion, enable_homology=True):
    """
    Parse the structure configuration data.

    Parameters
    ----------
    config : object
        The configuration data.
    time_explosion : float
        The time of the explosion.
    enable_homology : bool, optional
        Whether to enable homology (default is True).

    Returns
    -------
    electron_densities : object
        The parsed electron densities.
    temperature : object
        The parsed temperature.
    geometry : object
        The parsed geometry.
    density : object
        The parsed density.

    Raises
    ------
    NotImplementedError
        If the structure configuration type is not supported.

    Notes
    -----
    This function parses the structure configuration data and returns the parsed electron
    densities, temperature, geometry, and density. The structure configuration can be of
    type 'specific' or 'file'. If it is of type 'specific', the velocity and density are
    parsed from the configuration. If it is of type 'file', the velocity and density are
    read from a file. The parsed data is used to create a homologous radial 1D geometry object.
    """
    electron_densities = None
    temperature = None
    structure_config = config.model.structure
    if structure_config.type == "specific":
        velocity = quantity_linspace(
            structure_config.velocity.start,
            structure_config.velocity.stop,
            structure_config.velocity.num + 1,
        ).cgs
        density = parse_config_v1_density(config)

    elif structure_config.type == "file":
        if os.path.isabs(structure_config.filename):
            structure_config_fname = structure_config.filename
        else:
            structure_config_fname = os.path.join(
                config.config_dirname, structure_config.filename
            )

        (
            time_0,
            velocity,
            density_0,
            electron_densities,
            temperature,
        ) = read_density_file(structure_config_fname, structure_config.filetype)
        density_0 = density_0.insert(0, 0)

        density = calculate_density_after_time(
            density_0, time_0, time_explosion
        )

    else:
        raise NotImplementedError

    # Note: This is the number of shells *without* taking in mind the
    #       v boundaries.
    if len(density) == len(velocity):
        logger.warning(
            "Number of density points larger than number of shells. Assuming inner point irrelevant"
        )
        density = density[1:]
    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # v_inner
        velocity[1:],  # v_outer
        v_inner_boundary=structure_config.get("v_inner_boundary", None),
        v_outer_boundary=structure_config.get("v_outer_boundary", None),
        time_explosion=time_explosion,
    )
    return electron_densities, temperature, geometry, density
