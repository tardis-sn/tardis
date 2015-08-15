import numpy as np
import pandas as pd
from copy import deepcopy
from scipy import interpolate

from tardis.plasma.properties.base import (BasePlasmaProperty,
                                           ProcessingPlasmaProperty)

__all__ = ['PreviousElectronDensities', 'PreviousBetaSobolevs',
           'HeliumNLTE']

class PreviousIterationProperty(BasePlasmaProperty):
    def _set_output_value(self, output, value):
        setattr(self, output, value)

    def set_value(self, value):
        assert len(self.outputs) == 1
        self._set_output_value(self.outputs[0], value)

class PreviousElectronDensities(PreviousIterationProperty):
    outputs = ('previous_electron_densities',)

class PreviousBetaSobolevs(PreviousIterationProperty):
    outputs = ('previous_beta_sobolevs',)

class HeliumNLTE(ProcessingPlasmaProperty):
    outputs = ('helium_population',)

    def calculate(self, level_boltzmann_factor, electron_densities,
        ionization_data, beta_rad, g, g_electron, w, t_rad, t_electrons,
        delta, zeta_data, number_density):
        helium_population = deepcopy(level_boltzmann_factor.ix[2])
        # He I excited states
        he_one_population = level_boltzmann_factor.ix[2].ix[0] * \
            electron_densities * np.exp(
            ionization_data.ionization_energy.ix[2].ix[1] * beta_rad) * \
            (1 / (2 * g.ix[2].ix[1].ix[0] * g_electron * w) * (t_rad /
            t_electrons)**0.5)
        helium_population.ix[0].update(he_one_population)
        #He I metastable states
        helium_population.ix[0].ix[1] *= (1 / w)
        helium_population.ix[0].ix[2] *= (1 / w)
        #He I ground state
        helium_population.ix[0].ix[0] = 0.0
        #He II excited states
        he_two_population = level_boltzmann_factor.ix[2,1].mul(w * \
            (g.ix[2,1].ix[0]**(-1)) * (t_electrons / t_rad)**0.5)
        helium_population.ix[1].update(he_two_population)
        #He II ground state
        helium_population.ix[1].ix[0] = 1.0
        #He III states
        zeta = interpolate.interp1d(zeta_data.columns.values,
            zeta_data.ix[2,2].values)(t_rad)
        delta = delta.ix[2].ix[2]
        he_three_population = 2 * ((electron_densities)**(-1)) * g_electron * \
            (t_electrons / t_rad)**0.5 * np.exp(-1 *
            ionization_data.ionization_energy.ix[2,2] * beta_rad) * w * \
            ((delta * zeta) + w * (1 - zeta))
        helium_population.ix[0].ix[0].update(he_three_population)
        unnormalised = helium_population.sum()
        normalised = helium_population.mul(1.0 / unnormalised)
        helium_population.update(normalised)
        return helium_population

