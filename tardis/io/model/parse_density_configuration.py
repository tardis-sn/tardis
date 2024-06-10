from typing import Tuple

import numpy as np
from astropy import units as u

from tardis.io.configuration.config_reader import (
    Configuration,
    ConfigurationNameSpace,
)
from tardis.util.base import quantity_linspace


def parse_density_section(
    density_configuration: ConfigurationNameSpace,
    v_middle: u.Quantity,
    time_explosion: u.Quantity,
) -> Tuple[u.Quantity, u.Quantity]:
    """
    Parse the density section of the configuration file and produce a density at
    time_explosion.

    Parameters
    ----------
    density_configuration : tardis.io.config_reader.Configuration
    v_middle : astropy.Quantity
        middle of the velocity bins
    time_explosion : astropy.Quantity
        time of the explosion

    Returns
    -------
    density_0 : astropy.Quantity
        density at time_0
    time_0 : astropy.Quantity
        time of the density profile
    """
    if density_configuration.type == "branch85_w7":
        density_0 = calculate_power_law_density(
            v_middle,
            density_configuration.w7_v_0,
            density_configuration.w7_rho_0,
            -7,
        )
        time_0 = density_configuration.w7_time_0
    elif density_configuration.type == "uniform":
        density_0 = density_configuration.value.to("g cm^-3") * np.ones_like(
            v_middle.value
        )
        time_0 = density_configuration.get("time_0", time_explosion)
    elif density_configuration.type == "power_law":
        density_0 = calculate_power_law_density(
            v_middle,
            density_configuration.v_0,
            density_configuration.rho_0,
            density_configuration.exponent,
        )
        time_0 = density_configuration.get("time_0", time_explosion)
    elif density_configuration.type == "exponential":
        density_0 = calculate_exponential_density(
            v_middle, density_configuration.v_0, density_configuration.rho_0
        )
        time_0 = density_configuration.get("time_0", time_explosion)
    else:
        raise ValueError(
            f"Unrecognized density type '{density_configuration.type}'"
        )
    return density_0, time_0


def parse_config_v1_density(config: Configuration) -> u.Quantity:
    """
    Parse the configuration file and produce a density at
    time_explosion.

    Parameters
    ----------
    config : tardis.io.config_reader.Configuration

    Returns
    -------
    density: u.Quantity

    """
    velocity = quantity_linspace(
        config.model.structure.velocity.start,
        config.model.structure.velocity.stop,
        config.model.structure.velocity.num + 1,
    ).cgs

    adjusted_velocity = velocity.insert(0, 0)
    v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
    time_explosion = config.supernova.time_explosion.cgs
    d_conf = config.model.structure.density
    density_0, time_0 = parse_density_section(d_conf, v_middle, time_explosion)

    return calculate_density_after_time(density_0, time_0, time_explosion)


def parse_csvy_density(
    csvy_model_config: Configuration, time_explosion: u.Quantity
) -> u.Quantity:
    """
    Parse the density section of the csvy file and produce a density at
    time_explosion.

    Parameters
    ----------
    config : tardis.io.config_reader.Configuration
    csvy_model_config : tardis.io.config_reader.Configuration

    Returns
    -------
    density: u.Quantity

    """
    if hasattr(csvy_model_config, "velocity"):
        velocity = quantity_linspace(
            csvy_model_config.velocity.start,
            csvy_model_config.velocity.stop,
            csvy_model_config.velocity.num + 1,
        ).cgs
    else:
        velocity_field_index = [
            field.name for field in csvy_model_config.datatype.fields
        ].index("velocity")
        velocity_unit = u.Unit(
            csvy_model_config.datatype.fields[velocity_field_index].unit
        )
        velocity = csvy_model_config.velocity.values * velocity_unit

    adjusted_velocity = velocity.insert(0, 0)
    v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
    no_of_shells = len(adjusted_velocity) - 1

    if hasattr(csvy_model_config, "density"):
        density_0, time_0 = parse_density_section(
            csvy_model_config.density, v_middle, time_explosion
        )
    return calculate_density_after_time(density_0, time_0, time_explosion)


def calculate_power_law_density(
    velocities: u.Quantity,
    velocity_0: u.Quantity,
    rho_0: u.Quantity,
    exponent: float,
) -> u.Quantity:
    """

    This function computes a descret exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\left( \\frac{v}{v_0} \\right)^n`

    Parameters
    ----------
    velocities : astropy.Quantity
        Array like velocity profile
    velocity_0 : astropy.Quantity
        reference velocity
    rho_0 : astropy.Quantity
        reference density
    exponent : float
        exponent used in the powerlaw

    Returns
    -------
    densities : astropy.Quantity

    """
    densities = rho_0 * np.power((velocities / velocity_0), exponent)
    return densities


def calculate_exponential_density(
    velocities: u.Quantity, velocity_0: u.Quantity, rho_0: u.Quantity
) -> u.Quantity:
    """
    This function computes the exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\exp \\left( -\\frac{v}{v_0} \\right)`

    Parameters
    ----------
    velocities : astropy.Quantity
        Array like velocity profile
    velocity_0 : astropy.Quantity
        reference velocity
    rho_0 : astropy.Quantity
        reference density

    Returns
    -------
    densities : astropy.Quantity

    """
    densities = rho_0 * np.exp(-(velocities / velocity_0))
    return densities


def calculate_density_after_time(
    densities: u.Quantity, time_0: u.Quantity, time_explosion: u.Quantity
) -> u.Quantity:
    """
    scale the density from an initial time of the model to the
    time of the explosion by ^-3

    Parameters
    ----------
    densities : astropy.units.Quantity
        densities
    time_0 : astropy.units.Quantity
        time of the model
    time_explosion : astropy.units.Quantity
        time to be scaled to

    Returns
    -------
    scaled_density : astropy.units.Quantity
        in g / cm^3
    """
    return (densities * (time_explosion / time_0) ** -3).to(u.g / (u.cm**3))


def parse_density_csvy(csvy_model_config, csvy_model_data, time_explosion):
    """
    Parse the density data from a CSVY model.

    Parameters
    ----------
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    time_explosion : float
        The time of the explosion.

    Returns
    -------
    density : object
        The parsed density data.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the density data from a CSVY model. If the CSVY model configuration
    contains a 'density' attribute, it uses the 'parse_csvy_density' function to parse the
    density data. Otherwise, it calculates the density data using the 'calculate_density_after_time'
    function. The parsed density data is returned.
    """
    if hasattr(csvy_model_config, "density"):
        density = parse_csvy_density(csvy_model_config, time_explosion)
    else:
        time_0 = csvy_model_config.model_density_time_0
        density_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("density")
        density_unit = u.Unit(
            csvy_model_config.datatype.fields[density_field_index]["unit"]
        )
        density_0 = csvy_model_data["density"].values * density_unit
        # Removing as thee new architecture removes the 0th shell already
        # density_0 = density_0.to("g/cm^3")[1:]
        # density_0 = density_0.insert(0, 0)
        density = calculate_density_after_time(
            density_0, time_0, time_explosion
        )

    return density
