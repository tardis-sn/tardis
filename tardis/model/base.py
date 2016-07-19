import numpy as np
from astropy import constants, units as u

from tardis.util import quantity_linspace
from tardis.io.config_reader import calculate_power_law_density


class HomologousDensity(object):
    def __init__(self, density_0, time_0):
        self.density_0 = density_0
        self.time_0 = time_0

    def after_time(self, time_explosion):
        return (self.density_0 * (time_explosion / self.time_0) ** -3).cgs

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


class Radial1DModel(object):
    def __init__(self, velocity, homologous_density, abundance, t_inner,
                 time_explosion, t_radiative=None, dilution_factor=None,
                 v_boundary_inner=None, v_boundary_outer=None):
        self._v_boundary_inner = None
        self._v_boundary_outer = None
        self._velocity = None
        self.v_boundary_inner = v_boundary_inner
        self.v_boundary_outer = v_boundary_outer
        self.raw_velocity = velocity
        self.homologous_density = homologous_density
        self._abundance = abundance
        self.t_inner = t_inner
        self.time_explosion = time_explosion

        if t_radiative is None:
            lambda_wien_inner = constants.b_wien / self.t_inner
            self.t_radiative = constants.b_wien / (lambda_wien_inner * (
                1 + (self.v_middle - self.v_boundary_inner) / constants.c))
        else:
            self.t_radiative = t_radiative

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
    def radius(self):
        return self.time_explosion * self.velocity

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
    def velocity(self):
        if not self._velocity:
            self._velocity = self.raw_velocity[self.v_boundary_inner_index:
                                        self.v_boundary_outer_index].copy()
            self._velocity[0] = self.v_boundary_inner
            self._velocity[-1] = self.v_boundary_outer
        return self._velocity

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
        density = self.homologous_density.after_time(self.time_explosion)
        return density[self.v_boundary_inner_index
                       :self.v_boundary_outer_index][1:]

    @property
    def abundance(self):
        abundance = self._abundance.ix[:, self.v_boundary_inner_index:
                                     self.v_boundary_outer_index]
        abundance.columns = range(len(abundance.columns))
        return abundance

    @property
    def volume(self):
        return ((4. / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)).cgs

    @property
    def no_of_shells(self):
        return len(self.velocity) - 1

    @property
    def v_boundary_inner(self):
        if self._v_boundary_inner is None:
            return self.raw_velocity[0]
        return self._v_boundary_inner

    @v_boundary_inner.setter
    def v_boundary_inner(self, value):
        if value is not None:
            value = u.Quantity(value, self.v_boundary_inner.unit)
            if value > self.v_boundary_outer:
                raise ValueError('v_boundary_inner must not be higher than '
                                 'v_boundary_outer.')
            if value > self.raw_velocity[-1]:
                raise ValueError('v_boundary_inner is outside of '
                                 'the model range.')
            if value <= self.raw_velocity[0]:
                value = None
        self._v_boundary_inner = value
        # Invalidate the cached cut-down velocity array
        self._velocity = None

    @property
    def v_boundary_outer(self):
        if self._v_boundary_outer is None:
            return self.raw_velocity[-1]
        return self._v_boundary_outer

    @v_boundary_outer.setter
    def v_boundary_outer(self, value):
        if value is not None:
            value = u.Quantity(value, self.v_boundary_outer.unit)
            if value < self.v_boundary_inner:
                raise ValueError('v_boundary_outer must not be smaller than '
                                 'v_boundary_inner.')
            if value < self.raw_velocity[0]:
                raise ValueError('v_boundary_outer is outside of '
                                 'the model range.')
            if value >= self.raw_velocity[-1]:
                value = None
        self._v_boundary_outer = value
        # Invalidate the cached cut-down velocity array
        self._velocity = None

    @property
    def v_boundary_inner_index(self):
        if self.v_boundary_inner <= self.raw_velocity[0]:
            return None
        else:
            idx = max(0,
                      self.raw_velocity.searchsorted(self.v_boundary_inner) - 1)
            # check for zero volume of designated first cell
            if np.isclose(self.v_boundary_inner, self.raw_velocity[idx + 1],
                          atol=1e-8 * u.km / u.s) and (self.v_boundary_inner <=
                                                           self.raw_velocity[idx + 1]):
                idx += 1
            return idx

    @property
    def v_boundary_outer_index(self):
        if self.v_boundary_outer >= self.raw_velocity[-1]:
            return None
        return self.raw_velocity.searchsorted(self.v_boundary_outer) + 1

    @classmethod
    def from_config(cls, config):
        # TODO: Currently only supporting structure type specific
        t_inner = config.plasma.t_inner
        t_radiative = config.plasma.t_rads
        time_explosion = config.supernova.time_explosion
        abundance = config.abundances

        structure = config.model.structure
        if structure.type == 'specific':
            velocity = quantity_linspace(structure.velocity.start,
                                         structure.velocity.stop,
                                         structure.velocity.num + 1)
            homologous_density = HomologousDensity.from_config(config)
        else:
            raise NotImplementedError

        return cls(velocity=velocity,
                   homologous_density=homologous_density,
                   abundance=abundance,
                   t_inner=t_inner,
                   time_explosion=time_explosion,
                   t_radiative=t_radiative,
                   dilution_factor=None,
                   v_boundary_inner=None,
                   v_boundary_outer=None)
