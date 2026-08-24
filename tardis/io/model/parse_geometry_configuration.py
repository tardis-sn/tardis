from pathlib import Path

import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.readers.base import read_density_file
from tardis.model.geometry.radial1d import Radial1DGeometry
from tardis.model.geometry.radial1d_homologous import HomologousRadial1DGeometry
from tardis.util.base import quantity_linspace


def parse_structure_from_config(config: Configuration):
    """Parses the structure section from a config object

    Parameters
    ----------
    config
        The configuration to parse

    Returns
    -------
    density_time
        Time at which densities are valid.
    velocity
        Velocities.
    density
        Densities.
    electron_densities
        Electron densities.
    temperature
        Temperatures.

    Raises
    ------
    NotImplementedError
        For structure types that are not "specific" or "file"
    """
    density_time = None
    velocity = None
    density = None
    electron_densities = None
    temperature = None
    structure_config = config.model.structure
    if structure_config.type == "specific":
        velocity = quantity_linspace(
            structure_config.velocity.start,
            structure_config.velocity.stop,
            structure_config.velocity.num + 1,
        ).cgs

    elif structure_config.type == "file":
        if Path(structure_config.filename).is_absolute():
            structure_config_fname = structure_config.filename
        else:
            structure_config_fname = (
                Path(config.config_dirname) / structure_config.filename
            )

        (
            density_time,
            velocity,
            density,
            electron_densities,
            temperature,
        ) = read_density_file(structure_config_fname, structure_config.filetype)
    else:
        raise NotImplementedError

    return density_time, velocity, density, electron_densities, temperature


def parse_geometry_from_config(config: Configuration, time_explosion):
    """Parse the geometry data from a TARDIS config.

    Parameters
    ----------
    config
        Configuration object.
    time_explosion
        The time of the explosion.

    Returns
    -------
    HomologousRadial1DGeometry
        The parsed homologous geometry.
    """
    _, velocity, _, _, _ = parse_structure_from_config(config)

    return HomologousRadial1DGeometry(
        velocity[:-1],
        velocity[1:],
        v_inner_boundary=config.model.structure.get("v_inner_boundary", None),
        v_outer_boundary=config.model.structure.get("v_outer_boundary", None),
        time_explosion=time_explosion,
    )


def parse_nonhomologous_geometry_from_config(
    config: Configuration,
) -> Radial1DGeometry:
    """Parse nonhomologous geometry from a TARDIS configuration.

    Parameters
    ----------
    config : Configuration
        Configuration containing explicit velocity and radius boundaries.

    Returns
    -------
    Radial1DGeometry
        Geometry retaining independent radius and velocity boundaries.
    """
    _, velocity, _, _, _ = parse_structure_from_config(config)
    structure_config = config.model.structure
    radius = quantity_linspace(
        structure_config.radius.start,
        structure_config.radius.stop,
        structure_config.velocity.num + 1,
    ).to("cm")

    return Radial1DGeometry(
        r_inner=radius[:-1],
        r_outer=radius[1:],
        v_inner=velocity[:-1],
        v_outer=velocity[1:],
        r_inner_boundary=radius[0],
        r_outer_boundary=radius[-1],
        v_inner_boundary=structure_config.get("v_inner_boundary", None),
        v_outer_boundary=structure_config.get("v_outer_boundary", None),
    )


def parse_velocity_from_csvy(
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame | None,
) -> u.Quantity:
    """Parse shell-boundary velocities from a CSVY model.

    Parameters
    ----------
    csvy_model_config : Configuration
        Configuration data from the CSVY model.
    csvy_model_data : pandas.DataFrame or None
        Tabular data from the CSVY model.

    Returns
    -------
    astropy.units.Quantity
        Shell-boundary velocities in CGS units.
    """
    if hasattr(csvy_model_config, "velocity"):
        return quantity_linspace(
            csvy_model_config.velocity.start,
            csvy_model_config.velocity.stop,
            csvy_model_config.velocity.num + 1,
        ).cgs

    velocity_field_idx = [
        field["name"] for field in csvy_model_config.datatype.fields
    ].index("velocity")
    velocity_unit = u.Unit(
        csvy_model_config.datatype.fields[velocity_field_idx]["unit"]
    )
    velocity = csvy_model_data["velocity"].values * velocity_unit
    return velocity.to("cm/s")


def parse_homologous_geometry_from_csvy(
    config: Configuration,
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame | None,
    time_explosion: u.Quantity,
) -> HomologousRadial1DGeometry:
    """Parse homologous geometry from a CSVY model.

    Parameters
    ----------
    config : Configuration
        Main TARDIS configuration data.
    csvy_model_config : Configuration
        Configuration data from the CSVY model.
    csvy_model_data : pandas.DataFrame or None
        Tabular data from the CSVY model.
    time_explosion : astropy.units.Quantity
        Time since explosion.

    Returns
    -------
    HomologousRadial1DGeometry
        Homologous geometry constructed from the CSVY velocities.
    """
    model = getattr(config, "model", None)
    velocity = parse_velocity_from_csvy(csvy_model_config, csvy_model_data)

    return HomologousRadial1DGeometry(
        velocity[:-1],
        velocity[1:],
        v_inner_boundary=getattr(model, "v_inner_boundary", None),
        v_outer_boundary=getattr(model, "v_outer_boundary", None),
        time_explosion=time_explosion,
    )


def parse_nonhomologous_geometry_from_csvy(
    config: Configuration,
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame,
) -> Radial1DGeometry:
    """Parse nonhomologous geometry from a CSVY model.

    Parameters
    ----------
    config : Configuration
        Main TARDIS configuration data.
    csvy_model_config : Configuration
        Configuration data from the CSVY model.
    csvy_model_data : pandas.DataFrame
        Tabular radius and velocity data from the CSVY model.

    Returns
    -------
    Radial1DGeometry
        Nonhomologous geometry retaining the supplied radii and velocities.
    """
    model = getattr(config, "model", None)
    v_inner_boundary = getattr(model, "v_inner_boundary", None)
    v_outer_boundary = getattr(model, "v_outer_boundary", None)
    velocity = parse_velocity_from_csvy(csvy_model_config, csvy_model_data)

    radius_field_idx = [
        field["name"] for field in csvy_model_config.datatype.fields
    ].index("radius")
    radius_unit = u.Unit(
        csvy_model_config.datatype.fields[radius_field_idx]["unit"]
    )
    radius = (csvy_model_data["radius"].values * radius_unit).to("cm")

    return Radial1DGeometry(
        r_inner=radius[:-1],
        r_outer=radius[1:],
        v_inner=velocity[:-1],
        v_outer=velocity[1:],
        r_inner_boundary=radius[0],
        r_outer_boundary=radius[-1],
        v_inner_boundary=(
            v_inner_boundary if v_inner_boundary is not None else velocity[0]
        ),
        v_outer_boundary=(
            v_outer_boundary if v_outer_boundary is not None else velocity[-1]
        ),
    )
