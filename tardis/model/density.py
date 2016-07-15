"""These Density classes make use of descriptors, with the assumption that
~tardis.util.InstanceDescriptorMixin is used to achieve having instance
descriptors.

These classes are specifically designed to be attributes of a
model instance. When accessed from the model they dynamically calculate and
return the density array according to various values of the model, such as the
velocity, or the time_explosion.

Therefore, changing the velocity array, or the time_explosion in the model
will automagically change the density array too."""

import numpy as np


class Density(object):
    @staticmethod
    def after_time(density_0, time_0, time_explosion):
        if time_0 is None:
            time_0 = time_explosion
        else:
            time_0 = time_0
        return density_0 * (time_explosion / time_0) ** -3

    @classmethod
    def from_config(cls, config):
        d_conf = config.structure.density
        if d_conf.type == 'uniform':
            return UniformDensity(d_conf.value.to('g cm^-3'))
        elif d_conf.type == 'branch85_w7':
            return Branch85w7Density(d_conf.w7_time_0, d_conf.w7_rho_0,
                                     d_conf.w7_v_0)
        elif d_conf.type == 'power_law':
            return PowerLawDensity(d_conf.time_0, d_conf.rho_0,
                                   d_conf.v_0, d_conf.exponent)
        elif d_conf.type == 'exponential':
            return ExponentialDensity(d_conf.time_0, d_conf.rho_0, d_conf.v_0)


class StaticDensity(Density):
    """For density read from a file, doesn't change."""

    def __init__(self, density_0):
        self.density_0 = density_0

    def __get__(self, instance, owner):
        return self.density_0.cgs


class UniformDensity(Density):
    """One density value across the shells"""

    def __init__(self, density_0):
        self.density_0 = density_0

    def __get__(self, instance, owner):
        return np.ones(instance.no_of_shells) * self.density_0.cgs


class PowerLawDensity(Density):
    """For both power_law and branch85_w7 densities"""

    def __init__(self, time_0, rho_0, v_0, exponent):
        self.time_0 = time_0
        self.rho_0 = rho_0
        self.v_0 = v_0
        self.exponent = exponent

    def __get__(self, instance, owner):
        density_0 = self.calculate_power_law_density(instance.v_middle,
                                                     self.v_0, self.rho_0,
                                                     self.exponent)
        return (self.after_time(density_0, self.time_0,
                                instance.time_explosion)).cgs

    @staticmethod
    def calculate_power_law_density(velocities, velocity_0, rho_0, exponent):
        """
        This function computes a discrete exponential density profile.
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


class Branch85w7Density(PowerLawDensity):
    """"""

    def __init__(self, w7_time_0, w7_rho_0, w7_v_0):
        super(Branch85w7Density, self).__init__(w7_time_0, w7_rho_0, w7_v_0, -7)


class ExponentialDensity(Density):
    """For exponential density"""

    def __init__(self, time_0, rho_0, v_0):
        self.time_0 = time_0
        self.rho_0 = rho_0
        self.v_0 = v_0

    def __get__(self, instance, owner):
        density_0 = self.calculate_exponential_density(instance.v_middle,
                                                       self.v_0, self.rho_0)
        return (self.after_time(density_0, self.time_0,
                               instance.time_explosion)).cgs

    @staticmethod
    def calculate_exponential_density(velocities, v_0, rho0):
        """
        This function computes the exponential density profile.
        :math:`\\rho = \\rho_0 \\times \\exp \\left( -\\frac{v}{v_0} \\right)`

        Parameters
        ----------

        velocities : ~astropy.Quantity
            Array like velocity profile
        velocity_0 : ~astropy.Quantity
            reference velocity
        rho0 : ~astropy.Quantity
            reference density

        Returns
        -------

        densities : ~astropy.Quantity

        """
        densities = rho0 * np.exp(-(velocities / v_0))
        return densities
