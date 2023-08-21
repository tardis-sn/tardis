import logging
import os

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
)
from tardis.plasma.properties.ion_population import PhiSahaNebular

__all__ = [
    "HeliumNLTE",
    "HeliumNumericalNLTE",
]

logger = logging.getLogger(__name__)


class HeliumNLTE(ProcessingPlasmaProperty):
    outputs = ("helium_population",)

    @staticmethod
    def calculate(
        level_boltzmann_factor,
        ionization_data,
        beta_rad,
        g,
        g_electron,
        w,
        t_rad,
        t_electrons,
        delta,
        zeta_data,
        number_density,
        partition_function,
    ):
        """
        Updates all of the helium level populations according to the helium NLTE recomb approximation.
        """
        helium_population = level_boltzmann_factor.loc[2].copy()
        # He I excited states
        he_one_population = HeliumNLTE.calculate_helium_one(
            g_electron, beta_rad, ionization_data, level_boltzmann_factor, g, w
        )
        helium_population.loc[0].update(he_one_population)
        # He I ground state
        helium_population.loc[0, 0] = 0.0
        # He II excited states
        he_two_population = level_boltzmann_factor.loc[2, 1].mul(
            (float(g.loc[2, 1, 0]) ** (-1.0))
        )
        helium_population.loc[1].update(he_two_population)
        # He II ground state
        helium_population.loc[1, 0] = 1.0
        # He III states
        helium_population.loc[2, 0] = HeliumNLTE.calculate_helium_three(
            t_rad,
            w,
            zeta_data,
            t_electrons,
            delta,
            g_electron,
            beta_rad,
            ionization_data,
            g,
        )
        #        unnormalised = helium_population.sum()
        #        normalised = helium_population.mul(number_density.ix[2] / unnormalised)
        #        helium_population.update(normalised)
        return helium_population

    @staticmethod
    def calculate_helium_one(
        g_electron, beta_rad, ionization_data, level_boltzmann_factor, g, w
    ):
        """
        Calculates the He I level population values, in equilibrium with the He II ground state.
        """
        return (
            level_boltzmann_factor.loc[2, 0]
            * (1.0 / (2 * float(g.loc[2, 1, 0])))
            * (1 / g_electron)
            * (1 / (w**2.0))
            * np.exp(ionization_data.loc[2, 1] * beta_rad)
        )

    @staticmethod
    def calculate_helium_three(
        t_rad,
        w,
        zeta_data,
        t_electrons,
        delta,
        g_electron,
        beta_rad,
        ionization_data,
        g,
    ):
        """
        Calculates the He III level population values.
        """
        zeta = PhiSahaNebular.get_zeta_values(zeta_data, 2, t_rad)[1]
        he_three_population = (
            2
            * (float(g.loc[2, 2, 0]) / float(g.loc[2, 1, 0]))
            * g_electron
            * np.exp(-ionization_data.loc[2, 2] * beta_rad)
            * w
            * (delta.loc[2, 2] * zeta + w * (1.0 - zeta))
            * (t_electrons / t_rad) ** 0.5
        )
        return he_three_population


class HeliumNumericalNLTE(ProcessingPlasmaProperty):
    """
    IMPORTANT: This particular property requires a specific numerical NLTE
    solver and a specific atomic dataset (neither of which are distributed
    with Tardis) to work.
    """

    outputs = ("helium_population",)

    def __init__(self, plasma_parent, heating_rate_data_file):
        super(HeliumNumericalNLTE, self).__init__(plasma_parent)
        self._g_upper = None
        self._g_lower = None
        self.heating_rate_data = np.loadtxt(heating_rate_data_file, unpack=True)

    def calculate(
        self,
        ion_number_density,
        electron_densities,
        t_electrons,
        w,
        lines,
        j_blues,
        levels,
        level_boltzmann_factor,
        t_rad,
        zeta_data,
        g_electron,
        delta,
        partition_function,
        ionization_data,
        beta_rad,
        g,
        time_explosion,
    ):
        logger.info("Performing numerical NLTE He calculations.")
        if len(j_blues) == 0:
            return None
        # Outputting data required by SH module
        for zone, _ in enumerate(electron_densities):
            with open(
                f"He_NLTE_Files/shellconditions_{zone}.txt", "w"
            ) as output_file:
                output_file.write(ion_number_density.loc[2].sum()[zone])
                output_file.write(electron_densities[zone])
                output_file.write(t_electrons[zone])
                output_file.write(self.heating_rate_data[zone])
                output_file.write(w[zone])
                output_file.write(time_explosion)
                output_file.write(t_rad[zone])
                output_file.write(self.plasma_parent.v_inner[zone])
                output_file.write(self.plasma_parent.v_outer[zone])

        for zone, _ in enumerate(electron_densities):
            with open(
                f"He_NLTE_Files/abundances_{zone}.txt", "w"
            ) as output_file:
                for element in range(1, 31):
                    try:
                        number_density = (
                            ion_number_density[zone].loc[element].sum()
                        )
                    except:
                        number_density = 0.0
                        logger.debug(
                            f"Number Density could not be calculated. Setting Number Density to {number_density}"
                        )
                    output_file.write(number_density)

            helium_lines = lines[lines["atomic_number"] == 2]
            helium_lines = helium_lines[helium_lines["ion_number"] == 0]
        for zone, _ in enumerate(electron_densities):
            with open(
                f"He_NLTE_Files/discradfield_{zone}.txt", "w"
            ) as output_file:
                j_blues = pd.DataFrame(j_blues, index=lines.index)
                helium_j_blues = j_blues[zone].loc[helium_lines.index]
                for value in helium_lines.index:
                    if helium_lines.level_number_lower.loc[value] < 35:
                        output_file.write(
                            int(helium_lines.level_number_lower.loc[value] + 1),
                            int(helium_lines.level_number_upper.loc[value] + 1),
                            j_blues[zone].loc[value],
                        )
        # Running numerical simulations
        for zone, _ in enumerate(electron_densities):
            os.rename(
                f"He_NLTE_Files/abundances{zone}.txt",
                "He_NLTE_Files/abundances_current.txt",
            )
            os.rename(
                f"He_NLTE_Files/shellconditions{zone}.txt",
                "He_NLTE_Files/shellconditions_current.txt",
            )
            os.rename(
                f"He_NLTE_Files/discradfield{zone}.txt",
                "He_NLTE_Files/discradfield_current.txt",
            )
            os.system("nlte-solver-module/bin/nlte_solvertest >/dev/null")
            os.rename(
                "He_NLTE_Files/abundances_current.txt",
                f"He_NLTE_Files/abundances{zone}.txt",
            )
            os.rename(
                "He_NLTE_Files/shellconditions_current.txt",
                f"He_NLTE_Files/shellconditions{zone}.txt",
            )
            os.rename(
                "He_NLTE_Files/discradfield_current.txt",
                f"He_NLTE_Files/discradfield{zone}.txt",
            )
            os.rename("debug_occs.dat", f"He_NLTE_Files/occs{zone}.txt")
        # Reading in populations from files
        helium_population = level_boltzmann_factor.loc[2].copy()
        for zone, _ in enumerate(electron_densities):
            with open(
                f"He_NLTE_Files/discradfield{zone}.txt", "r"
            ) as read_file:
                for level in range(0, 35):
                    level_population = read_file.readline()
                    level_population = float(level_population)
                    helium_population[zone].loc[0, level] = level_population
                helium_population[zone].loc[1, 0] = float(read_file.readline())
        # Performing He LTE level populations (upper two energy levels,
        # He II excited states, He III)
        he_one_population = HeliumNLTE.calculate_helium_one(
            g_electron,
            beta_rad,
            partition_function,
            ionization_data,
            level_boltzmann_factor,
            electron_densities,
            g,
            w,
            t_rad,
            t_electrons,
        )
        helium_population.loc[0, 35].update(he_one_population.loc[35])
        helium_population.loc[0, 36].update(he_one_population.loc[36])

        he_two_population = level_boltzmann_factor.loc[2, 1, 1:].mul(
            (float(g.loc[2, 1, 0]) ** (-1)) * helium_population.loc[s1, 0]
        )
        helium_population.loc[1, 1:].update(he_two_population)

        helium_population.loc[2, 0] = HeliumNLTE.calculate_helium_three(
            t_rad,
            w,
            zeta_data,
            t_electrons,
            delta,
            g_electron,
            beta_rad,
            partition_function,
            ionization_data,
            electron_densities,
        )
        unnormalised = helium_population.sum()
        normalised = helium_population.mul(
            ion_number_density.loc[2].sum() / unnormalised
        )
        helium_population.update(normalised)
        return helium_population
