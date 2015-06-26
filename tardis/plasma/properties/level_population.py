import logging

import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelPopulation', 'LevelNumberDensity',
'LevelBoltzmannFactorNLTE']


class LevelPopulation(ProcessingPlasmaProperty):

        name = 'level_population_fraction'
        latexformula = (r'$N_{i, j, k} = \frac{g_{i, j, k} '
                         r'e^{-\beta_\textrm{rad} E_{i, j, k}}}{Z_{i, j}}$')

        @staticmethod
        def calculate(levels, partition_function, level_boltzmann_factor):

            partition_function_broadcast = partition_function.ix[
                levels.index.droplevel(2)].values

            return level_boltzmann_factor / partition_function_broadcast

class LevelNumberDensity(ProcessingPlasmaProperty):
    name = 'level_number_density'

    @staticmethod
    def calculate(level_population_fraction, ion_number_density):
        ion_number_density_broadcast = ion_number_density.ix[
            level_population_fraction.index.droplevel(2)].values
        return level_population_fraction * ion_number_density_broadcast

class LevelBoltzmannFactorNLTE(ProcessingPlasmaProperty):

    name = 'level_boltzmann_factor_nlte'

    @staticmethod
    def calculate(nlte_input, beta_sobolev, j_blues, levels,
        nlte_data, t_electron, electron_densities,
        general_level_boltzmann_factor):
        """
        Calculating the NLTE level populations for specific ions

        """

        if nlte_input.get('coronal_approximation', False):
            beta_sobolev_nlte = np.ones_like(beta_sobolev)
            j_blues_nlte = np.zeros_like(j_blues)
            logger.info('using coronal approximation = setting beta_sobolevs to 1 AND j_blues to 0')
        else:
            beta_sobolev_nlte = beta_sobolev
            j_blues_nlte = j_blues.values

        if nlte_input.get('classical_nebular', False):
            logger.info('using Classical Nebular = setting beta_sobolevs to 1')
            beta_sobolev_nlte = np.ones_like(beta_sobolev)

        for species in nlte_input.species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = len(levels.ix[species])

            lnl = nlte_data.lines_level_number_lower[species]
            lnu = nlte_data.lines_level_number_upper[species]

            lines_index = nlte_data.lines_idx[species]
            A_uls = nlte_data.A_uls[species]
            B_uls = nlte_data.B_uls[species]
            B_lus = nlte_data.B_lus[species]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            r_ul_matrix = np.zeros((number_of_levels,
                number_of_levels, len(t_electron)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + \
                B_uls[np.newaxis].T * j_blues_nlte[lines_index]
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolev_nlte[lines_index]

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
                j_blues_nlte[lines_index] * beta_sobolev_nlte[lines_index]

            collision_matrix = nlte_data.get_collision_matrix(
                species, t_electron) * electron_density.values

            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

            rates_matrix[0, :, :] = 1.0

            level_boltzmann_factor = general_level_boltzmann_factor

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(t_electron)):
                nlte_level_boltzmann_factor = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                level_boltzmann_factor[i].ix[species] = \
                    nlte_level_boltzmann_factor
            return level_boltzmann_factor