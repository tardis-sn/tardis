import numpy as np
import pandas as pd
from scipy import interpolate

from tardis.plasma.properties.base import (BasePlasmaProperty,
                                           ProcessingPlasmaProperty)
from tardis.plasma.properties import PhiSahaNebular, PhiSahaLTE

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
        delta, zeta_data, number_density, partition_function):
        helium_population = level_boltzmann_factor.ix[2].copy()
        # He I excited states
        he_one_population = self.calculate_helium_one(g_electron, beta_rad,
            partition_function, ionization_data, level_boltzmann_factor,
            electron_densities, g, w, t_rad, t_electrons)
        helium_population.ix[0].update(he_one_population)
        #He I metastable states
        helium_population.ix[0].ix[1] *= (1 / w)
        helium_population.ix[0].ix[2] *= (1 / w)
        #He I ground state
        helium_population.ix[0].ix[0] = 0.0
        #He II excited states
        he_two_population = level_boltzmann_factor.ix[2,1].mul(
            (g.ix[2,1].ix[0]**(-1)) * (t_electrons / t_rad)**0.5)
        helium_population.ix[1].update(he_two_population)
        #He II ground state
        helium_population.ix[1].ix[0] = 1.0
        #He III states
        helium_population.ix[2].ix[0] = self.calculate_helium_three(t_rad, w,
            zeta_data, t_electrons, delta, g_electron, beta_rad,
            partition_function, ionization_data, electron_densities)
        unnormalised = helium_population.sum()
        normalised = helium_population.mul(number_density.ix[2] / unnormalised)
        helium_population.update(normalised)
        return helium_population

    def calculate_helium_one(self, g_electron, beta_rad, partition_function,
            ionization_data, level_boltzmann_factor, electron_densities, g,
            w, t_rad, t_electron):
        (partition_function_index, ionization_data_index, partition_function,
            ionization_data) = self.filter_with_helium_index(2, 1,
            partition_function, ionization_data)
        phis = (1 / PhiSahaLTE.calculate(g_electron, beta_rad,
            partition_function, ionization_data)) * electron_densities * \
            (1.0/g.ix[2].ix[1].ix[0]) * (1/w) * (t_rad/t_electron)**(0.5)
        return level_boltzmann_factor.ix[2].ix[0].mul(
            pd.DataFrame(phis.ix[2].ix[1].values)[0].transpose())

    def calculate_helium_three(self, t_rad, w, zeta_data, t_electrons, delta,
        g_electron, beta_rad, partition_function, ionization_data,
        electron_densities):
        (partition_function_index, ionization_data_index, partition_function,
            ionization_data) = self.filter_with_helium_index(2, 2,
            partition_function, ionization_data)
        zeta_data = pd.DataFrame(zeta_data.ix[2].ix[2].values,
            columns=ionization_data_index, index=zeta_data.columns).transpose()
        delta = pd.DataFrame(delta.ix[2].ix[2].values,
            columns=ionization_data_index, index=delta.columns).transpose()
        phis = PhiSahaNebular.calculate(t_rad, w,
            zeta_data, t_electrons, delta, g_electron,
            beta_rad, partition_function, ionization_data)
        return (phis * (partition_function.ix[2].ix[1] /
            partition_function.ix[2].ix[2]) * (1 /
            electron_densities)).ix[2].ix[2]

    @staticmethod
    def filter_with_helium_index(atomic_number, ion_number, partition_function,
        ionization_data):
        partition_function_index = pd.MultiIndex.from_tuples([(atomic_number,
            ion_number-1), (atomic_number, ion_number)],
            names=['atomic_number', 'ion_number'])
        ionization_data_index = pd.MultiIndex.from_tuples([(atomic_number,
            ion_number)],
            names=['atomic_number', 'ion_number'])
        partition_function = pd.DataFrame(
            partition_function.ix[atomic_number].ix[
            ion_number-1:ion_number].values,
            index=partition_function_index, columns=partition_function.columns)
        ionization_data = pd.DataFrame(
            ionization_data.ix[atomic_number].ix[ion_number][
            'ionization_energy'], index=ionization_data_index,
            columns=['ionization_energy'])
        return partition_function_index, ionization_data_index,\
               partition_function, ionization_data