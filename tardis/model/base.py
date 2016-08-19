import os
import logging
import numpy as np
import pandas as pd
from astropy import constants, units as u

from tardis.util import quantity_linspace, element_symbol2atomic_number
from tardis.io.model_reader import read_density_file, read_abundances_file
from tardis.io.util import to_hdf
from density import HomologousDensity

logger = logging.getLogger(__name__)


class Radial1DModel(object):
    def __init__(self, velocity, homologous_density, abundance, time_explosion,
                 t_inner, luminosity_requested=None, t_radiative=None,
                 dilution_factor=None, v_boundary_inner=None,
                 v_boundary_outer=None):
        self._v_boundary_inner = None
        self._v_boundary_outer = None
        self._velocity = None
        self.raw_velocity = velocity
        self.v_boundary_inner = v_boundary_inner
        self.v_boundary_outer = v_boundary_outer
        self.homologous_density = homologous_density
        self._abundance = abundance
        self.time_explosion = time_explosion
        if t_inner is None:
            if luminosity_requested is not None:
                self.t_inner = ((luminosity_requested /
                                (4 * np.pi * self.r_inner[0] ** 2 *
                                 constants.sigma_sb)) ** .25).to('K')
            else:
                raise ValueError('Both t_inner and luminosity_requested cannot '
                                 'be None.')
        else:
            self.t_inner = t_inner

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
        start = self.v_boundary_inner_index
        stop = self.v_boundary_outer_index
        if stop is not None:
            # abundance has one element less than velocity
            # ix stop index is inclusive
            stop -= 2
        abundance = self._abundance.ix[:, start:stop]
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

    def to_hdf(self, path_or_buf, path=''):
        """
        Store the model to an HDF structure.

        Parameters
        ----------
        path_or_buf
            Path or buffer to the HDF store
        path : str
            Path inside the HDF store to store the model

        Returns
        -------
        None

        """
        model_path = os.path.join(path, 'model')
        properties = ['t_inner', 'w', 't_radiative', 'v_inner', 'v_outer']
        to_hdf(path_or_buf, model_path, {name: getattr(self, name) for name
                                         in properties})

    @classmethod
    def from_config(cls, config):
        time_explosion = config.supernova.time_explosion.cgs

        structure = config.model.structure
        if structure.type == 'specific':
            velocity = quantity_linspace(structure.velocity.start,
                                         structure.velocity.stop,
                                         structure.velocity.num + 1).cgs
            homologous_density = HomologousDensity.from_config(config)
        elif structure.type == 'file':
            if os.path.isabs(structure.filename):
                structure_fname = structure.filename
            else:
                structure_fname = os.path.join(config.config_dirname,
                                               structure.filename)

            time_0, velocity, density_0 = read_density_file(
                structure_fname, structure.filetype)
            homologous_density = HomologousDensity(density_0, time_0)
        else:
            raise NotImplementedError
        no_of_shells = len(velocity) - 1

        if config.plasma.initial_t_rad > 0 * u.K:
            t_radiative = np.ones(no_of_shells) * config.plasma.initial_t_rad
        else:
            t_radiative = None

        if config.plasma.initial_t_inner < 0.0 * u.K:
            luminosity_requested = config.supernova.luminosity_requested
            t_inner = None
        else:
            luminosity_requested = None
            t_inner = config.plasma.initial_t_inner

        abundances_section = config.model.abundances
        if abundances_section.type == 'uniform':
            abundance = pd.DataFrame(columns=np.arange(no_of_shells),
                                     index=pd.Index(np.arange(1, 120),
                                                    name='atomic_number'),
                                     dtype=np.float64)

            for element_symbol_string in abundances_section:
                if element_symbol_string == 'type':
                    continue
                z = element_symbol2atomic_number(element_symbol_string)
                abundance.ix[z] = float(abundances_section[element_symbol_string])

        elif abundances_section.type == 'file':
            if os.path.isabs(abundances_section.filename):
                abundances_fname = abundances_section.filename
            else:
                abundances_fname = os.path.join(config.config_dirname,
                                                abundances_section.filename)

            index, abundance = read_abundances_file(abundances_fname,
                                                    abundances_section.filetype)

        abundance = abundance.replace(np.nan, 0.0)
        abundance = abundance[abundance.sum(axis=1) > 0]

        norm_factor = abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning("Abundances have not been normalized to 1."
                           " - normalizing")
            abundance /= norm_factor


        return cls(velocity=velocity,
                   homologous_density=homologous_density,
                   abundance=abundance,
                   time_explosion=time_explosion,
                   t_radiative=t_radiative,
                   t_inner=t_inner,
                   luminosity_requested=luminosity_requested,
                   dilution_factor=None,
                   v_boundary_inner=structure.get('v_inner_boundary', None),
                   v_boundary_outer=structure.get('v_outer_boundary', None))
