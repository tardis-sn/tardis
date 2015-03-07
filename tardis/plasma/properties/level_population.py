import logging

import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['LevelPopulationLTE', 'LevelNumberDensity']


class LevelPopulationLTE(ProcessingPlasmaProperty):

        name = 'level_population'
        latex_formula = (r'N_{i, j, k} = \frac{g_{i, j, k} '
                         r'e^{-\beta_\textrm{rad} E_{i, j, k}}{Z_{i, j}}')

        @staticmethod
        def calculate(levels, partition_function, level_boltzmann_factor):

            partition_function_broadcast = partition_function.ix[
                levels.index.droplevel(2)].values

            return level_boltzmann_factor / partition_function_broadcast

class LevelNumberDensity(ProcessingPlasmaProperty):
    name = 'level_number_density'

    @staticmethod
    def calculate(level_population, ion_number_density):
        ion_number_density_broadcast = ion_number_density.ix[
            level_population.index.droplevel(2)].values
        return level_population * ion_number_density_broadcast


class LevelPopulationDiluteLTE(ProcessingPlasmaProperty):
        """
        Calculate the level populations and putting them in the column 'number-density' of the self.levels table.
        :math:`N` denotes the ion number density calculated with `calculate_ionization_balance`, i is the atomic number,
        j is the ion number and k is the level number. For non-metastable levels we add the dilution factor (W) to the calculation.

        .. math::

            N_{i, j, k}(\\textrm{metastable}) &= \\frac{g_k}{Z_{i, j}}\\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} E_k} \\\\
            N_{i, j, k}(\\textrm{not metastable}) &= W\\frac{g_k}{Z_{i, j}}\\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} E_k} \\\\


        This function updates the 'number_density' column on the levels table (or adds it if non-existing)
        """

        @staticmethod
        def calculate(levels, partition_function,
                                        level_boltzmann_factor,
                                        ion_number_density, w):
            partition_function_broadcast = partition_function.ix[
                levels.index.droplevel(2)].values

            ion_number_density_broadcast = ion_number_density.ix[levels.index.droplevel(2)].values

            level_population = level_boltzmann_factor / partition_function_broadcast
            level_population[~levels.metastable] *= np.min([w, np.ones_like(w)],axis=0)


class LevelPopulationNLTE(ProcessingPlasmaProperty):
    @staticmethod
    def calculate(self):
        """
        Calculating the NLTE level populations for specific ions

        """

        if not hasattr(self, 'beta_sobolevs'):
            self.beta_sobolevs = np.zeros_like(self.tau_sobolevs.values)

        macro_atom.calculate_beta_sobolev(self.tau_sobolevs.values.ravel(order='F'),
                                          self.beta_sobolevs.ravel(order='F'))
        self.beta_sobolevs_precalculated = True

        if self.nlte_config.get('coronal_approximation', False):
            beta_sobolevs = np.ones_like(self.beta_sobolevs)
            j_blues = np.zeros_like(self.j_blues)
            logger.info('using coronal approximation = setting beta_sobolevs to 1 AND j_blues to 0')
        else:
            beta_sobolevs = self.beta_sobolevs
            j_blues = self.j_blues.values

        if self.nlte_config.get('classical_nebular', False):
            logger.info('using Classical Nebular = setting beta_sobolevs to 1')
            beta_sobolevs = np.ones_like(self.beta_sobolevs)

        for species in self.nlte_config.species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = self.atom_data.levels.energy.ix[species].count()

            level_populations = self.level_populations.ix[species].values
            lnl = self.atom_data.nlte_data.lines_level_number_lower[species]
            lnu = self.atom_data.nlte_data.lines_level_number_upper[species]

            lines_index = self.atom_data.nlte_data.lines_idx[species]
            A_uls = self.atom_data.nlte_data.A_uls[species]
            B_uls = self.atom_data.nlte_data.B_uls[species]
            B_lus = self.atom_data.nlte_data.B_lus[species]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            r_ul_matrix = np.zeros((number_of_levels, number_of_levels, len(self.t_rads)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2, len(self.t_rads)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + B_uls[np.newaxis].T * j_blues[lines_index]
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs[lines_index]

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2, len(self.t_rads)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * j_blues[lines_index] * beta_sobolevs[lines_index]

            collision_matrix = self.atom_data.nlte_data.get_collision_matrix(species, self.t_electrons) * \
                               self.electron_densities.values


            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

            rates_matrix[0, :, :] = 1.0

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(self.t_rads)):
                relative_level_populations = np.linalg.solve(rates_matrix[:, :, i], x)
                self.level_populations[i].ix[species] = relative_level_populations * self.ion_populations[i].ix[species]