from pathlib import Path

import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.readers.base import read_density_file
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.model.geometry.radial1d_nonhomologous import NonhomologousRadial1DGeometry
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
        density = density.insert(0, 0)
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
        The parsed geometry.
    """
    (
        density_time,
        velocity,
        density,
        electron_densities,
        temperature,
    ) = parse_structure_from_config(config)

    return HomologousRadial1DGeometry(
        velocity[:-1],  # v_inner
        velocity[1:],  # v_outer
        v_inner_boundary=config.model.structure.get("v_inner_boundary", None),
        v_outer_boundary=config.model.structure.get("v_outer_boundary", None),
        time_explosion=time_explosion,
    )


def parse_geometry_from_csvy(
    config: Configuration,
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame | None,
    time_explosion: u.Quantity,
) -> HomologousRadial1DGeometry | NonhomologousRadial1DGeometry:
    """Parse the geometry data from a CSVY model.

    Parameters
    ----------
    config
        The configuration data.
    csvy_model_config
        The configuration data of the CSVY model.
    csvy_model_data
        The data of the CSVY model.
    time_explosion
        The time of the explosion.

    Returns
    -------
    HomologousRadial1DGeometry
        The parsed geometry.

    Notes
    -----
    This function parses the geometry data from a CSVY model. It extracts the velocity
    information from the CSVY model configuration or data. The parsed velocity data is
    used to create a homologous radial 1D geometry object, which is returned.
    """
    if hasattr(config, "model"):
        if hasattr(config.model, "v_inner_boundary"):
            v_inner_boundary = config.model.v_inner_boundary
        else:
            v_inner_boundary = None

        if hasattr(config.model, "v_outer_boundary"):
            v_outer_boundary = config.model.v_outer_boundary
        else:
            v_outer_boundary = None
    else:
        v_inner_boundary = None
        v_outer_boundary = None

    if hasattr(csvy_model_config, "velocity"):
        velocity = quantity_linspace(
            csvy_model_config.velocity.start,
            csvy_model_config.velocity.stop,
            csvy_model_config.velocity.num + 1,
        ).cgs
    else:
        velocity_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("velocity")
        velocity_unit = u.Unit(
            csvy_model_config.datatype.fields[velocity_field_index]["unit"]
        )
        velocity = csvy_model_data["velocity"].values * velocity_unit
        velocity = velocity.to("cm/s")

    if csvy_model_data is not None and "radius" in csvy_model_data.columns:
        radius_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("radius")
        radius_unit = u.Unit(
            csvy_model_config.datatype.fields[radius_field_index]["unit"]
        )
        radius = (csvy_model_data["radius"].values * radius_unit).to("cm")
        return NonhomologousRadial1DGeometry(
            r_inner=radius[:-1],
            r_outer=radius[1:],
            v_inner=velocity[:-1],
            v_outer=velocity[1:],
            r_inner_boundary=radius[0],
            r_outer_boundary=radius[-1],
            v_inner_boundary=v_inner_boundary if v_inner_boundary is not None else velocity[0],
            v_outer_boundary=v_outer_boundary if v_outer_boundary is not None else velocity[-1],
        )
    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # v_inner
        velocity[1:],  # v_outer
        v_inner_boundary=v_inner_boundary,
        v_outer_boundary=v_outer_boundary,
        time_explosion=time_explosion,
    )
    return geometry
