import logging

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaConfigError

logger = logging.getLogger(__name__)

__all__ = ['LevelBoltzmannFactorLTE', 'LevelBoltzmannFactorDiluteLTE',
           'LevelBoltzmannFactorNoNLTE', 'LevelBoltzmannFactorNLTE',
           'PartitionFunction']

class LevelBoltzmannFactorLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    general_level_boltzmann_factor : Pandas DataFrame, dtype float
                             Level population proportionality values. Indexed
                             by atomic number, ion number, level number.
                             Columns corresponding to zones. Does not consider
                             NLTE.
    """
    outputs = ('general_level_boltzmann_factor',)
    latex_name = ('bf_{i,j,k}',)
    latex_formula = ('g_{i,j,k}e^{\\dfrac{-\\epsilon_{i,j,k}}{k_{\
        \\textrm{B}}T_{\\textrm{rad}}}}',)

    @staticmethod
    def calculate(excitation_energy, g, beta_rad, levels):
        exponential = np.exp(np.outer(excitation_energy.values, -beta_rad))
        level_boltzmann_factor_array = (g.values[np.newaxis].T *
                                        exponential)
        level_boltzmann_factor = pd.DataFrame(level_boltzmann_factor_array,
                                              index=levels,
                                              columns=np.arange(len(beta_rad)),
                                              dtype=np.float64)
        return level_boltzmann_factor

class LevelBoltzmannFactorDiluteLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    general_level_boltzmann_factor : Pandas DataFrame, dtype float
                             Level population proportionality values. Indexed
                             by atomic number, ion number, level number.
                             Columns corresponding to zones. Dilute radiation
                             field means non-metastable level values are
                             multiplied by an additional factor W. Does not
                             consider NLTE.
    """
    outputs = ('general_level_boltzmann_factor',)
    latex_name = ('bf_{i,j,k}',)
    latex_formula = ('Wg_{i,j,k}e^{\\dfrac{-\\epsilon_{i,j,k}}{k_{\
        \\textrm{B}}T_{\\textrm{rad}}}}',)

    def calculate(self, levels, g, excitation_energy, beta_rad, w,
        metastability):
        level_boltzmann_factor = LevelBoltzmannFactorLTE.calculate(
            excitation_energy, g, beta_rad, levels)
        level_boltzmann_factor[~metastability] *= w
        return level_boltzmann_factor

class LevelBoltzmannFactorNoNLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    level_boltzmann_factor : Pandas DataFrame, dtype float
                             Returns general_level_boltzmann_factor as this
                             property is included if NLTE is not used.
    """
    outputs = ('level_boltzmann_factor',)

    @staticmethod
    def calculate(general_level_boltzmann_factor):
        return general_level_boltzmann_factor

class LevelBoltzmannFactorNLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    level_boltzmann_factor : Pandas DataFrame, dtype float
                             Returns general_level_boltzmann_factor but
                             updated for those species treated in NLTE.
    """
    outputs = ('level_boltzmann_factor',)

    def calculate(self):
        pass

    def __init__(self, plasma_parent, classical_nebular=False,
        coronal_approximation=False):
        """
        Selects appropriate 'calculate' function based on NLTE config
        options selected.
        """
        super(LevelBoltzmannFactorNLTE, self).__init__(plasma_parent)
        if classical_nebular == True and coronal_approximation == False:
            self.calculate = self._calculate_classical_nebular
        elif coronal_approximation == True and classical_nebular == False:
            self.calculate = self._calculate_coronal_approximation
        elif coronal_approximation == True and classical_nebular == True:
            raise PlasmaConfigError('Both coronal approximation and classical'
                                    'nebular specified in the config.')
        else:
            self.calculate = self._calculate_general
        self._update_inputs()

    def _main_nlte_calculation(self, atomic_data, nlte_data,
        t_electrons, j_blues, beta_sobolevs, general_level_boltzmann_factor,
        previous_electron_densities):
        """
        The core of the NLTE calculation, used with all possible config.
        options.
        """
        for species in self.plasma_parent.nlte_species:
            j_blues = j_blues.values
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
                len(t_electrons)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2,
                len(t_electrons)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + \
                B_uls[np.newaxis].T * j_blues[lines_index]
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs[lines_index]
            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2,
                len(t_electrons)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
                j_blues[lines_index] * beta_sobolevs[lines_index]
            if atomic_data.has_collision_data:
                if previous_electron_densities is None:
                    collision_matrix = r_ul_matrix.copy()
                    collision_matrix.fill(0.0)
                else:
                    collision_matrix = nlte_data.get_collision_matrix(species,
                        t_electrons) * previous_electron_densities.values
            else:
                collision_matrix = r_ul_matrix.copy()
                collision_matrix.fill(0.0)
            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix
            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)
            rates_matrix[0, :, :] = 1.0
            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(t_electrons)):
                level_boltzmann_factor = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                general_level_boltzmann_factor[i].ix[species] = \
                    level_boltzmann_factor
        return general_level_boltzmann_factor

    def _calculate_classical_nebular(self, t_electrons, lines, atomic_data,
        nlte_data, general_level_boltzmann_factor, j_blues,
        lte_j_blues, previous_electron_densities):
        """
        Performs NLTE calculations using the classical nebular treatment.
        All beta sobolev values taken as 1.
        """
        beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        if len(j_blues)==0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(j_blues, index=lines.index, columns =
                range(len(t_electrons)))
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data, nlte_data, t_electrons, j_blues,
            beta_sobolevs, general_level_boltzmann_factor,
            previous_electron_densities)
        return general_level_boltzmann_factor

    def _calculate_coronal_approximation(self, t_electrons, lines, atomic_data,
        nlte_data, general_level_boltzmann_factor,
        previous_electron_densities):
        """
        Performs NLTE calculations using the coronal approximation.
        All beta sobolev values taken as 1 and j_blues taken as 0.
        """
        beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        j_blues = np.zeros((len(lines), len(t_electrons)))
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data, nlte_data, t_electrons, j_blues,
            beta_sobolevs, general_level_boltzmann_factor,
            previous_electron_densities)
        return general_level_boltzmann_factor

    def _calculate_general(self, t_electrons, lines, atomic_data, nlte_data,
        general_level_boltzmann_factor, j_blues,
        previous_beta_sobolev, lte_j_blues, previous_electron_densities):
        """
        Full NLTE calculation without approximations.
        """
        if previous_beta_sobolev is None:
            beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        else:
            beta_sobolevs = previous_beta_sobolev
        if len(j_blues)==0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(j_blues, index=lines.index, columns =
                range(len(t_electrons)))
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data, nlte_data, t_electrons, j_blues,
            beta_sobolevs, general_level_boltzmann_factor,
            previous_electron_densities)
        return general_level_boltzmann_factor

class PartitionFunction(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    partition_function : Pandas DataFrame, dtype float
                         Indexed by atomic number, ion number.
                         Columns are zones.
    """
    outputs = ('partition_function',)
    latex_name = ('Z_{i,j}',)
    latex_formula = ('\\sum_{k}bf_{i,j,k}',)

    def calculate(self, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=['atomic_number', 'ion_number']).sum()
