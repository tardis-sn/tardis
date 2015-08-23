import pandas as pd
import numpy as np
import os

from tardis.plasma.properties.base import (BasePlasmaProperty,
                                           ProcessingPlasmaProperty)
from tardis.plasma.properties import PhiSahaNebular, PhiSahaLTE

__all__ = ['PreviousElectronDensities', 'PreviousBetaSobolevs',
           'HeliumNLTE', 'HeliumNumericalNLTE']

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
            (g.ix[2,1].ix[0]**(-1)))
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

    @staticmethod
    def calculate_helium_one(g_electron, beta_rad, partition_function,
            ionization_data, level_boltzmann_factor, electron_densities, g,
            w, t_rad, t_electron):
        (partition_function_index, ionization_data_index, partition_function,
            ionization_data) = HeliumNLTE.filter_with_helium_index(2, 1,
            partition_function, ionization_data)
        phis = (1 / PhiSahaLTE.calculate(g_electron, beta_rad,
            partition_function, ionization_data)) * electron_densities * \
            (1.0/g.ix[2,1,0]) * (1/w) * (t_rad/t_electron)**(0.5)
        return level_boltzmann_factor.ix[2].ix[0].mul(
            pd.DataFrame(phis.ix[2].ix[1].values)[0].transpose())

    @staticmethod
    def calculate_helium_three(t_rad, w, zeta_data, t_electrons, delta,
        g_electron, beta_rad, partition_function, ionization_data,
        electron_densities):
        (partition_function_index, ionization_data_index, partition_function,
            ionization_data) = HeliumNLTE.filter_with_helium_index(2, 2,
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

class HeliumNumericalNLTE(ProcessingPlasmaProperty):
    outputs = ('helium_population',)
    '''
    IMPORTANT: This particular property requires a specific numerical NLTE
    solver and a specific atomic dataset (neither of which are distributed
    with Tardis) to work.
    '''
    def calculate(self, ion_number_density, electron_densities, t_electrons, w,
        lines, j_blues, levels, level_boltzmann_factor, t_rad,
        zeta_data, g_electron, delta, partition_function, ionization_data,
        beta_rad, g):
        print 'Performing numerical NLTE calculations.'
        if len(j_blues)==0:
            return None
        heating_rate_data = np.loadtxt(
            self.plasma_parent.heating_rate_data_file, unpack=True)
        #Outputting data required by SH module
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/shellconditions_%d.txt' %zone, 'w') as \
                    output_file:
                print>>output_file, ion_number_density.ix[2].sum()[zone]
                print>>output_file, electron_densities[zone]
                print>>output_file, t_electrons[zone]
                print>>output_file, heating_rate_data[zone]
                print>>output_file, w[zone]
                print>>output_file, self.plasma_parent.time_explosion
                print>>output_file, t_rad[zone]
                print>>output_file, self.plasma_parent.v_inner[zone]
                print>>output_file, self.plasma_parent.v_outer[zone]

        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/abundances_%d.txt' %zone, 'w') as \
                    output_file:
                for element in range(1,31):
                    try:
                        number_density = ion_number_density[zone].ix[
                            element].sum()
                    except:
                        number_density = 0.0
                    print>>output_file, number_density

            helium_lines = lines[lines['atomic_number']==2]
            helium_lines = helium_lines[helium_lines['ion_number']==0]
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/discradfield_%d.txt' %zone, 'w') as \
                    output_file:
                j_blues = pd.DataFrame(j_blues, index=lines.index)
                helium_j_blues = j_blues[zone].ix[helium_lines.index]
                for value in helium_lines.index:
                    if (helium_lines.level_number_lower.ix[value]<35):
                        print>>output_file, \
                            int(helium_lines.level_number_lower.ix[value]+1), \
                            int(helium_lines.level_number_upper.ix[value]+1), \
                            j_blues[zone].ix[value]
        #Running numerical simulations
        for zone, _ in enumerate(electron_densities):
            os.rename('He_NLTE_Files/abundances_%d.txt' %zone,
                      'He_NLTE_Files/abundances_current.txt')
            os.rename('He_NLTE_Files/shellconditions_%d.txt' %zone,
                      'He_NLTE_Files/shellconditions_current.txt')
            os.rename('He_NLTE_Files/discradfield_%d.txt' %zone,
                      'He_NLTE_Files/discradfield_current.txt')
            os.system("nlte-solver-module/bin/nlte_solvertest >/dev/null")
            os.rename('He_NLTE_Files/abundances_current.txt',
                      'He_NLTE_Files/abundances_%d.txt' %zone)
            os.rename('He_NLTE_Files/shellconditions_current.txt',
                      'He_NLTE_Files/shellconditions_%d.txt' %zone)
            os.rename('He_NLTE_Files/discradfield_current.txt',
                      'He_NLTE_Files/discradfield_%d.txt' %zone)
            os.rename('debug_occs.dat', 'He_NLTE_Files/occs_%d.txt' %zone)
        #Reading in populations from files
        helium_population = level_boltzmann_factor.ix[2].copy()
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/discradfield_%d.txt' %zone, 'r') as \
                    read_file:
                for level in range(0, 35):
                    level_population = read_file.readline()
                    level_population = float(level_population)
                    helium_population[zone].ix[0][level] = level_population
                helium_population[zone].ix[1].ix[0] = float(
                    read_file.readline())
        #Performing He LTE level populations (upper two energy levels,
        #He II excited states, He III)
        he_one_population = HeliumNLTE.calculate_helium_one(g_electron,
            beta_rad, partition_function, ionization_data,
            level_boltzmann_factor, electron_densities, g, w, t_rad,
            t_electrons)
        helium_population.ix[0].ix[35].update(he_one_population.ix[35])
        helium_population.ix[0].ix[36].update(he_one_population.ix[36])

        he_two_population = level_boltzmann_factor.ix[2].ix[1].ix[1:].mul(
            (g.ix[2,1,0]**(-1)) * helium_population.ix[s1,0])
        helium_population.ix[1].ix[1:].update(he_two_population)

        helium_population.ix[2].ix[0] = HeliumNLTE.calculate_helium_three(
            t_rad, w, zeta_data, t_electrons, delta, g_electron, beta_rad,
            partition_function, ionization_data, electron_densities)
        unnormalised = helium_population.sum()
        normalised = helium_population.mul(ion_number_density.ix[2].sum()
            / unnormalised)
        helium_population.update(normalised)
        return helium_population