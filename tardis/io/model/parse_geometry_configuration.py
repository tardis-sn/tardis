from astropy import units as u

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.util.base import quantity_linspace


def parse_csvy_geometry(
    config, csvy_model_config, csvy_model_data, time_explosion
):
    """
    Parse the geometry data from a CSVY model.

    Parameters
    ----------
    config : object
        The configuration data.
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    time_explosion : float
        The time of the explosion.

    Returns
    -------
    geometry : object
        The parsed geometry.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the geometry data from a CSVY model. It extracts the velocity
    information from the CSVY model configuration or data. The parsed velocity data is
    used to create a homologous radial 1D geometry object, which is returned.
    """
    if hasattr(config, "model"):
        if hasattr(config.model, "v_inner_boundary"):
            v_boundary_inner = config.model.v_inner_boundary
        else:
            v_boundary_inner = None

        if hasattr(config.model, "v_outer_boundary"):
            v_boundary_outer = config.model.v_outer_boundary
        else:
            v_boundary_outer = None
    else:
        v_boundary_inner = None
        v_boundary_outer = None

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

    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # v_inner
        velocity[1:],  # v_outer
        v_inner_boundary=v_boundary_inner,
        v_outer_boundary=v_boundary_outer,
        time_explosion=time_explosion,
    )
    return geometry
