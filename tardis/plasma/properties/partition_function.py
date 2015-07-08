import logging

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelBoltzmannFactorLTE', 'LevelBoltzmannFactorDiluteLTE',
           'LevelBoltzmannFactorNoNLTE', 'LevelBoltzmannFactorNLTECoronal',
           'LevelBoltzmannFactorNLTEGeneral', 'PartitionFunction']

class LevelBoltzmannFactorLTE(ProcessingPlasmaProperty):
    """
    Calculate the level population Boltzmann factor

    .. math:
        {latex_formula}

    """

    outputs = ('general_level_boltzmann_factor',)
    latex_formula = r'$g_{i, j, k} e^{E_{i, j, k} \times \beta_\textrm{rad}}$'

    @staticmethod
    def calculate(levels, beta_rad):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        return level_boltzmann_factor

class LevelBoltzmannFactorDiluteLTE(ProcessingPlasmaProperty):

    outputs = ('general_level_boltzmann_factor',)

    @staticmethod
    def calculate(levels, beta_rad, w):
        exponential = np.exp(np.outer(levels.energy.values, -beta_rad))
        level_boltzmann_factor_array = (levels.g.values[np.newaxis].T *
                                        exponential)

        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels.index,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        level_boltzmann_factor[~levels.metastable] *= w
        return level_boltzmann_factor

class LevelBoltzmannFactorNoNLTE(ProcessingPlasmaProperty):

    outputs = ('level_boltzmann_factor',)

    @staticmethod
    def calculate(general_level_boltzmann_factor):
        return general_level_boltzmann_factor

class LevelBoltzmannFactorNLTECoronal(ProcessingPlasmaProperty):

    outputs = ('level_boltzmann_factor',)

    @staticmethod
    def calculate(t_electron, lines, atomic_data, nlte_data,
            general_level_boltzmann_factor, nlte_species):
        """
        Calculating the NLTE level populations for specific ions

        """

        beta_sobolevs = np.ones((len(lines), len(t_electron)))
        j_blues = np.zeros((len(lines), len(t_electron)))

        for species in nlte_species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = atomic_data.levels.energy.ix[species].count()

            lnl = nlte_data.lines_level_number_lower[species]
            lnu = nlte_data.lines_level_number_upper[species]

            lines_index = nlte_data.lines_idx[species]
            A_uls = nlte_data.A_uls[species]
            B_uls = nlte_data.B_uls[species]
            B_lus = nlte_data.B_lus[species]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            r_ul_matrix = np.zeros((number_of_levels, number_of_levels,
                len(t_electron)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + \
                B_uls[np.newaxis].T * j_blues[lines_index]
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs[lines_index]

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
                j_blues[lines_index] * beta_sobolevs[lines_index]

#            collision_matrix = self.atom_data.nlte_data.get_collision_matrix(species, self.t_electrons) * \
#                               self.electron_densities.values
            collision_matrix = r_ul_matrix.copy()
            collision_matrix.fill(0.0)

            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

            rates_matrix[0, :, :] = 1.0

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(t_electron)):
                level_boltzmann_factor = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                general_level_boltzmann_factor[i].ix[species] = \
                    level_boltzmann_factor
        return general_level_boltzmann_factor

class LevelBoltzmannFactorNLTEClassicalNebular(ProcessingPlasmaProperty):
    pass

class LevelBoltzmannFactorNLTEGeneral(ProcessingPlasmaProperty):
    outputs = ('level_boltzmann_factor',)

    @staticmethod
    def calculate(t_electron, lines, atomic_data, nlte_data,
            general_level_boltzmann_factor, nlte_species, j_blues,
            previous_beta_sobolevs, lte_j_blues):
        """
        Calculating the NLTE level populations for specific ions

        """
        if previous_beta_sobolevs is None:
            beta_sobolevs = np.ones((len(lines), len(t_electron)))
        else:
            beta_sobolevs = previous_beta_sobolevs

        if len(j_blues)==0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(j_blues, index=lines.index, columns =
                range(len(t_electron)))

        for species in nlte_species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = atomic_data.levels.energy.ix[species].count()

            lnl = nlte_data.lines_level_number_lower[species]
            lnu = nlte_data.lines_level_number_upper[species]

            lines_index = nlte_data.lines_idx[species]
            A_uls = nlte_data.A_uls[species]
            B_uls = nlte_data.B_uls[species]
            B_lus = nlte_data.B_lus[species]

            j_blues_index = lines.index[lines_index]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            r_ul_matrix = np.zeros((number_of_levels, number_of_levels,
                len(t_electron)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + \
                B_uls[np.newaxis].T * j_blues.ix[j_blues_index]
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs[lines_index]

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
                j_blues.ix[j_blues_index] * beta_sobolevs[lines_index]

#            collision_matrix = self.atom_data.nlte_data.get_collision_matrix(species, self.t_electrons) * \
#                               self.electron_densities.values
            collision_matrix = r_ul_matrix.copy()
            collision_matrix.fill(0.0)

            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

            rates_matrix[0, :, :] = 1.0

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(t_electron)):
                level_boltzmann_factor = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                general_level_boltzmann_factor[i].ix[species] = \
                    level_boltzmann_factor
        return general_level_boltzmann_factor

class PartitionFunction(ProcessingPlasmaProperty):
    outputs = ('partition_function',)
    latex_outputs = '$Z_{i, j}$'

    latex_formula = (r'$Z_{i, j} = \sum_{k=1}^n g_{i, j, k} '
                     r'e^{E_{i, j, k} \times \beta_\textrm{rad}}$')

    @staticmethod
    def calculate(levels, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=['atomic_number', 'ion_number']).sum()
