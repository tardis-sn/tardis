import numpy as np


class HomologousDensity(object):
    def __init__(self, density_0, time_0):
        self.density_0 = density_0
        self.time_0 = time_0

    def after_time(self, time_explosion):
        return self.density_0 * (time_explosion / self.time_0) ** -3


class Radial1DModel(object):
    def __init__(self, velocity, homologous_density, abundance, t_inner,
                 t_radiative, v_boundary_inner, v_boundary_outer,
                 time_explosion, dilution_factor=None):
        self.velocity = velocity
        self.homologous_density = homologous_density
        self.abundance = abundance
        self.t_inner = t_inner
        self.t_radiative = t_radiative
        self.v_boundary_inner = v_boundary_inner
        self.v_boundary_outer = v_boundary_outer
        self.time_explosion = time_explosion
        if dilution_factor is None:
            self.dilution_factor = 0.5 * (1 - np.sqrt(
                1 - (self.r_inner[0] ** 2 / self.r_middle ** 2).to(1).value))
        else:
            self.dilution_factor = dilution_factor

    @property
    def w(self):
        return self.dilution_factor

    @w.setter
    def w(self, value):
        self.dilution_factor = value

    @property
    def t_rad(self):
        return self.t_radiative

    @t_rad.setter
    def t_rad(self, value):
        self.t_radiative = value

    @property
    def r_inner(self):
        return self.time_explosion * self.v_inner

    @property
    def r_outer(self):
        return self.time_explosion * self.v_outer

    @property
    def r_middle(self):
        return 0.5 * self.r_inner + 0.5 * self.r_outer

    @property
    def v_inner(self):
        return self.velocity[:-1]

    @property
    def v_outer(self):
        return self.velocity[1:]

    @property
    def v_middle(self):
        return 0.5 * self.v_inner + 0.5 * self.v_outer

    @property
    def density(self):
        return self.homologous_density.after_time(self.time_explosion)

    @property
    def volume(self):
        return ((4. / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)).cgs

    @property
    def no_of_shells(self):
        return len(self.velocity) - 1

    @classmethod
    def from_config(cls, config):
        pass
