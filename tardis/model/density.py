import numpy as np

from tardis.util import quantity_linspace


class HomologousDensity(object):
    def __init__(self, density_0, time_0):
        self.density_0 = density_0
        self.time_0 = time_0

    def after_time(self, time_explosion):
        return calculate_density_after_time(self.density_0, self.time_0,
                                            time_explosion).cgs

    @classmethod
    def from_config(cls, config):
        d_conf = config.structure.density
        if d_conf.type == 'branch85_w7':
            # This is temporary, until the old model is removed.
            velocity = quantity_linspace(config.structure.velocity.start,
                                         config.structure.velocity.stop,
                                         config.structure.velocity.num + 1)
            adjusted_velocity = velocity.insert(0, 0)
            v_middle = (adjusted_velocity[1:] * 0.5 +
                        adjusted_velocity[:-1] * 0.5)
            density_0 = calculate_power_law_density(v_middle, d_conf.w7_v_0,
                                                    d_conf.w7_rho_0, -7)
            return cls(density_0, d_conf.w7_time_0)
        else:
            raise NotImplementedError


def calculate_power_law_density(velocities, velocity_0, rho_0, exponent):
    """

    This function computes a descret exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\left( \\frac{v}{v_0} \\right)^n`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho_0 : ~astropy.Quantity
        reference density
    exponent : ~float
        exponent used in the powerlaw

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho_0 * np.power((velocities / velocity_0), exponent)
    return densities


def calculate_exponential_density(velocities, velocity_0, rho_0):
    """
    This function computes the exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\exp \\left( -\\frac{v}{v_0} \\right)`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho_0 : ~astropy.Quantity
        reference density

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho_0 * np.exp(-(velocities / velocity_0))
    return densities


def calculate_density_after_time(densities, time_0, time_explosion):
    """
    scale the density from an initial time of the model to the
    time of the explosion by ^-3

    Parameters:
    -----------

    densities: ~astropy.units.Quantity
        densities

    time_0: ~astropy.units.Quantity
        time of the model

    time_explosion: ~astropy.units.Quantity
        time to be scaled to

    Returns:
    --------

    scaled_density
    """

    return densities * (time_explosion / time_0) ** -3
