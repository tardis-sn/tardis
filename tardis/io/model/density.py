import numpy as np
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.util.base import quantity_linspace


def parse_density_section(
    density_configuration: Configuration,
    v_middle: u.Quantity,
    time_explosion: u.Quantity,
):
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
            f"Unrecognized density type " f"'{density_configuration.type}'"
        )
    return density_0, time_0


def parse_config_v1_density(config: Configuration) -> u.Quantity:
    """
    Create a new HomologousDensity instance from a Configuration object.

    Parameters
    ----------
    config : tardis.io.config_reader.Configuration

    Returns
    -------
    HomologousDensity

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


def parse_csvy_density(csvy_model_config, time_explosion: u.Quantity):
    """
    Create a new HomologousDensity instance from a base
    Configuration object and a csvy model Configuration object.

    Parameters
    ----------
    config : tardis.io.config_reader.Configuration
    csvy_model_config : tardis.io.config_reader.Configuration

    Returns
    -------
    HomologousDensity

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


def calculate_power_law_density(velocities, velocity_0, rho_0, exponent):
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


def calculate_exponential_density(velocities, velocity_0, rho_0):
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
