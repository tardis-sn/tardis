import numpy as np
from astropy import units as u
import logging

from tardis.util.base import quantity_linspace
from tardis.io.util import HDFWriterMixin, config_iteratation

logger = logging.getLogger(__name__)


class HomologousDensity(HDFWriterMixin):
    """A class that holds an initial density and time

    Parameters
    ----------
    density_0 : astropy.units.Quantity
    time_0 : astropy.units.Quantity

    """

    hdf_properties = ["density_0", "time_0"]

    def __init__(self, density_0, time_0):
        self.density_0 = density_0
        self.time_0 = time_0

    def calculate_density_at_time_of_simulation(self, time_explosion):
        """
        Scale `density_0` from an `time_0` to the `time_explosion` by ^-3

        Parameters
        ----------
        time_explosion : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        return calculate_density_after_time(
            self.density_0, self.time_0, time_explosion
        ).cgs

    @classmethod
    def from_csvy(cls, config, csvy_model_config):
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
        time_explosion = config.supernova.time_explosion.cgs

        if hasattr(csvy_model_config, "density"):
            d_conf = csvy_model_config.density
            density_type = d_conf.type
            if density_type == "branch85_w7":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.w7_v_0, d_conf.w7_rho_0, -7
                )
                time_0 = d_conf.w7_time_0
            elif density_type == "uniform":
                density_0 = d_conf.value.to("g cm^-3") * np.ones(no_of_shells)
                time_0 = d_conf.get("time_0", time_explosion)
            elif density_type == "power_law":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.v_0, d_conf.rho_0, d_conf.exponent
                )
                time_0 = d_conf.get("time_0", time_explosion)
            elif density_type == "exponential":
                density_0 = calculate_exponential_density(
                    v_middle, d_conf.v_0, d_conf.rho_0
                )
                time_0 = d_conf.get("time_0", time_explosion)
            else:
                raise ValueError(
                    f"Unrecognized density type " f"'{d_conf.type}'"
                )
        return cls(density_0, time_0)

    @classmethod
    def from_config(cls, config):
        """
        Create a new HomologousDensity instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        Returns
        -------
        HomologousDensity

        """
        d_conf = config.model.structure.density
        velocity = quantity_linspace(
            config.model.structure.velocity.start,
            config.model.structure.velocity.stop,
            config.model.structure.velocity.num + 1,
        ).cgs

        adjusted_velocity = velocity.insert(0, 0)
        v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        no_of_shells = len(adjusted_velocity) - 1
        time_explosion = config.supernova.time_explosion.cgs

        if d_conf.type == "branch85_w7":
            density_0 = calculate_power_law_density(
                v_middle, d_conf.w7_v_0, d_conf.w7_rho_0, -7
            )
            time_0 = d_conf.w7_time_0
        elif d_conf.type == "uniform":
            density_0 = d_conf.value.to("g cm^-3") * np.ones(no_of_shells)
            time_0 = d_conf.get("time_0", time_explosion)
        elif d_conf.type == "power_law":
            density_0 = calculate_power_law_density(
                v_middle, d_conf.v_0, d_conf.rho_0, d_conf.exponent
            )
            time_0 = d_conf.get("time_0", time_explosion)
        elif d_conf.type == "exponential":
            density_0 = calculate_exponential_density(
                v_middle, d_conf.v_0, d_conf.rho_0
            )
            time_0 = d_conf.get("time_0", time_explosion)
        else:
            raise ValueError(f"Unrecognized density type " f"'{d_conf.type}'")
        log_string = config_iteratation(d_conf)
        logger.debug(f"Density Config:\n\t  {log_string}")
        return cls(density_0, time_0)


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


def calculate_density_after_time(densities, time_0, time_explosion):
    """
    scale the density from an initial time of the model to the
    time of the explosion by ^-3

    Parameters
    -----------
    densities : astropy.units.Quantity
        densities
    time_0 : astropy.units.Quantity
        time of the model
    time_explosion : astropy.units.Quantity
        time to be scaled to

    Returns
    --------
    scaled_density
    """

    return densities * (time_explosion / time_0) ** -3
