import logging
import os

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import (PreviousIterationProperty,
                                           ProcessingPlasmaProperty)
from tardis.plasma.properties import PhiSahaNebular, PhiSahaLTE

__all__ = ['PreviousElectronDensities', 'PreviousBetaSobolev',
           'HeliumNLTE', 'HeliumNumericalNLTE']

logger = logging.getLogger(__name__)

class PreviousElectronDensities(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_electron_densities : The values for the electron densities converged upon in the previous iteration.
    """
    outputs = ('previous_electron_densities',)

    def set_initial_value(self, kwargs):
        initial_value = np.ones(len(kwargs['abundance'].columns))*1000000.0
        self._set_initial_value(initial_value)

class PreviousBetaSobolev(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_beta_sobolev : The beta sobolev values converged upon in the previous iteration.
    """
    outputs = ('previous_beta_sobolev',)

    def set_initial_value(self, kwargs):
        try:
            lines = len(kwargs['atomic_data'].lines)
        except:
            lines = len(kwargs['atomic_data']._lines)
        initial_value = np.ones((lines,
            len(kwargs['abundance'].columns)))
        self._set_initial_value(initial_value)

class HeliumNLTE(ProcessingPlasmaProperty):
    outputs = ('helium_population',)

    def calculate(self, level_boltzmann_factor, electron_densities,
        ionization_data, beta_rad, g, g_electron, w, t_rad, t_electrons,
        delta, zeta_data, number_density, partition_function):
        """
        Updates all of the helium level populations according to the helium NLTE recomb approximation.
        """
        helium_population = level_boltzmann_factor.ix[2].copy()
        # He I excited states
        he_one_population = self.calculate_helium_one(g_electron, beta_rad,
            ionization_data, level_boltzmann_factor, electron_densities, g, w)
        helium_population.ix[0].update(he_one_population)
        #He I metastable states
        helium_population.ix[0,1] *= (1 / w)
        helium_population.ix[0,2] *= (1 / w)
        #He I ground state
        helium_population.ix[0,0] = 0.0
        #He II excited states
        he_two_population = level_boltzmann_factor.ix[2,1].mul(
            (g.ix[2,1].ix[0]**(-1)))
        helium_population.ix[1].update(he_two_population)
        #He II ground state
        helium_population.ix[1,0] = 1.0
        #He III states
        helium_population.ix[2,0] = self.calculate_helium_three(t_rad, w,
            zeta_data, t_electrons, delta, g_electron, beta_rad,
            ionization_data, electron_densities, g)
        unnormalised = helium_population.sum()
        normalised = helium_population.mul(number_density.ix[2] / unnormalised)
        helium_population.update(normalised)
        return helium_population

    @staticmethod
    def calculate_helium_one(g_electron, beta_rad, ionization_data,
        level_boltzmann_factor, electron_densities, g, w):
        """
        Calculates the He I level population values, in equilibrium with the He II ground state.
        """
        return level_boltzmann_factor.ix[2,0].mul(
            g.ix[2,0], axis=0) * (1./(2*g.ix[2,1,0])) * \
            (1/g_electron) * (1/(w**2)) * np.exp(
            ionization_data.ionization_energy.ix[2,1] * beta_rad) * \
            electron_densities

    @staticmethod
    def calculate_helium_three(t_rad, w, zeta_data, t_electrons, delta,
        g_electron, beta_rad, ionization_data, electron_densities, g):
        """
        Calculates the He III level population values.
        """
        zeta = PhiSahaNebular.get_zeta_values(zeta_data, 2, t_rad)[1]
        he_three_population = (2 / electron_densities) * \
            (float(g.ix[2,2,0])/g.ix[2,1,0]) * g_electron * \
            np.exp(-ionization_data.ionization_energy.ix[2,2] * beta_rad) \
            * w * (delta.ix[2,2] * zeta + w * (1. - zeta)) * \
            (t_electrons / t_rad) ** 0.5

class HeliumNumericalNLTE(ProcessingPlasmaProperty):
    '''
    IMPORTANT: This particular property requires a specific numerical NLTE
    solver and a specific atomic dataset (neither of which are distributed
    with Tardis) to work.
    '''
    outputs = ('helium_population',)
    def calculate(self, ion_number_density, electron_densities, t_electrons, w,
        lines, j_blues, levels, level_boltzmann_factor, t_rad,
        zeta_data, g_electron, delta, partition_function, ionization_data,
        beta_rad, g):
        logger.info('Performing numerical NLTE He calculations.')
        if len(j_blues)==0:
            return None
        heating_rate_data = np.loadtxt(
            self.plasma_parent.heating_rate_data_file, unpack=True)
        #Outputting data required by SH module
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/shellconditions_{}.txt'.format(zone),
                'w') as output_file:
                output_file.write(ion_number_density.ix[2].sum()[zone])
                output_file.write(electron_densities[zone])
                output_file.write(t_electrons[zone])
                output_file.write(heating_rate_data[zone])
                output_file.write(w[zone])
                output_file.write(self.plasma_parent.time_explosion)
                output_file.write(t_rad[zone])
                output_file.write(self.plasma_parent.v_inner[zone])
                output_file.write(self.plasma_parent.v_outer[zone])

        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/abundances_{}.txt'.format(zone), 'w') as \
                    output_file:
                for element in range(1,31):
                    try:
                        number_density = ion_number_density[zone].ix[
                            element].sum()
                    except:
                        number_density = 0.0
                    output_file.write(number_density)

            helium_lines = lines[lines['atomic_number']==2]
            helium_lines = helium_lines[helium_lines['ion_number']==0]
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/discradfield_{}.txt'.format(zone), 'w') \
                    as output_file:
                j_blues = pd.DataFrame(j_blues, index=lines.index)
                helium_j_blues = j_blues[zone].ix[helium_lines.index]
                for value in helium_lines.index:
                    if (helium_lines.level_number_lower.ix[value]<35):
                        output_file.write(
                            int(helium_lines.level_number_lower.ix[value]+1),
                            int(helium_lines.level_number_upper.ix[value]+1),
                            j_blues[zone].ix[value])
        #Running numerical simulations
        for zone, _ in enumerate(electron_densities):
            os.rename('He_NLTE_Files/abundances{}.txt'.format(zone),
                      'He_NLTE_Files/abundances_current.txt')
            os.rename('He_NLTE_Files/shellconditions{}.txt'.format(zone),
                      'He_NLTE_Files/shellconditions_current.txt')
            os.rename('He_NLTE_Files/discradfield{}.txt'.format(zone),
                      'He_NLTE_Files/discradfield_current.txt')
            os.system("nlte-solver-module/bin/nlte_solvertest >/dev/null")
            os.rename('He_NLTE_Files/abundances_current.txt',
                      'He_NLTE_Files/abundances{}.txt'.format(zone))
            os.rename('He_NLTE_Files/shellconditions_current.txt',
                      'He_NLTE_Files/shellconditions{}.txt'.format(zone))
            os.rename('He_NLTE_Files/discradfield_current.txt',
                      'He_NLTE_Files/discradfield{}.txt'.format(zone))
            os.rename('debug_occs.dat', 'He_NLTE_Files/occs{}.txt'.format(zone))
        #Reading in populations from files
        helium_population = level_boltzmann_factor.ix[2].copy()
        for zone, _ in enumerate(electron_densities):
            with open('He_NLTE_Files/discradfield{}.txt'.format(zone), 'r') as \
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
