import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData


class NLTEExcitationSolver(ProcessingPlasmaProperty):
    def calculate(
        self,
        atomic_data,
        j_blues,
        beta_sobolev,
        nlte_excitation_species,
        gamma,
        alpha_sp,
        alpha_stim,
        coll_ion_coeff,
        coll_recomb_coeff,
        coll_exc_coeff,
        coll_deexc_coeff,
        electron_density,
    ):
        nlte_data = NLTEExcitationData(
            atomic_data.lines, nlte_excitation_species
        )
        rate_matrix_blocks = dict.fromkeys(nlte_excitation_species)
        for species in nlte_excitation_species:
            number_of_levels = atomic_data.levels.energy.loc[species].count()
            (
                lines_index,
                r_ul_index,
                r_ul_matrix,
                r_lu_index,
                r_lu_matrix,
            ) = self.prepare_r_uls_r_lus(
                number_of_levels, 1, j_blues, species, nlte_data
            )
            bound_rate_matrix = self.prepare_bound_bound_rate_matrix(
                number_of_levels,
                lines_index,
                r_ul_index,
                r_ul_matrix,
                r_lu_index,
                r_lu_matrix,
                beta_sobolev,
            )
            

    @staticmethod
    def prepare_bound_bound_rate_matrix(
        number_of_levels,
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
        beta_sobolev,
    ):
        """Calculates a matrix with bound-bound rates for NLTE excitation treatment.

        Parameters
        ----------
        number_of_levels : int
            Number of levels for the specified species.
        lines_index : numpy.array
            Index of lines in nlte_data.
        r_ul_index : numpy.array
            Index used for r_ul matrix
        r_ul_matrix : numpy.array
            Matrix with the rates(upper to lower transition) of bound-bound interaction(DOES NOT INCLUDE THE BETA SOBOLEVS)
            (number_of_levels, number_of_levels, number_of_shells)
        r_lu_index : numpy.array
            Index used for r_lu matrix
        r_lu_matrix : numpy.array
        r_lu_matrix : numpy.array
            Matrix with the rates(lower to upper transition) of bound-bound interaction(DOES NOT INCLUDE THE BETA SOBOLEVS)
            (number_of_levels, number_of_levels, number_of_shells)
        beta_sobolev : pandas.DataFrame
            Beta Sobolev factors.
        Returns
        -------
        numpy.array (number of levels, number of levels)
            Matrix with excitation-deexcitation rates(should be added to NLTE rate matrix for excitation treatment).
            NOTE: only works with ONE ion number treated in NLTE excitation AT ONCE.
        """
        number_of_shells = beta_sobolev.shape[1]
        try:
            beta_sobolev_filtered = beta_sobolev.iloc[lines_index]
        except AttributeError:
            beta_sobolev_filtered = beta_sobolev
        r_ul_matrix_reshaped = r_ul_matrix.reshape(
            (number_of_levels**2, number_of_shells)
        )
        r_lu_matrix_reshaped = r_lu_matrix.reshape(
            (number_of_levels**2, number_of_shells)
        )
        r_ul_matrix_reshaped[r_ul_index] *= beta_sobolev_filtered
        r_lu_matrix_reshaped[r_lu_index] *= beta_sobolev_filtered
        rates_matrix_bound_bound = r_ul_matrix + r_lu_matrix
        for i in range(number_of_levels):
            rates_matrix_bound_bound[i, i] = -rates_matrix_bound_bound[
                :, i
            ].sum(axis=0)
        return rates_matrix_bound_bound

    @staticmethod
    def prepare_r_uls_r_lus(
        number_of_levels,
        number_of_shells,
        j_blues,
        excitation_species,
        nlte_data,
    ):
        """Calculates part of rate matrices for bound bound interactions.

        Parameters
        ----------
        number_of_levels : int
            Number of levels for the NLTE excitation species.
        number_of_shells : int
            Number of shells.
        j_blues : pandas.DataFrame, dtype float
            Mean intensities in the blue wings of the line transitions.
        excitation_species : tuple
            Species treated in NLTE excitation.
        nlte_data : NLTEExcitationData
            Data relevant to NLTE excitation species.

        Returns
        -------
        lines_index : numpy.array
            Index of lines in nlte_data.
        number_of_levels : int
            Number of levels for the specified species.
        r_ul_index : numpy.array
            Index used for r_ul matrix
        r_ul_matrix_reshaped : numpy.array
            Matrix with the rates(upper to lower transition) of bound-bound interaction(DOES NOT INCLUDE THE BETA SOBOLEVS)
        r_lu_index : numpy.array
            Index used for r_lu matrix
        r_lu_matrix_reshaped : numpy.array
            Matrix with the rates(lower to upper transition) of bound-bound interaction(DOES NOT INCLUDE THE BETA SOBOLEVS)
        """
        # number_of_levels = atomic_data_levels.energy.loc[
        #     excitation_species
        # ].count() do this in the solver
        lnl = nlte_data.lines_level_number_lower[excitation_species]
        lnu = nlte_data.lines_level_number_upper[excitation_species]
        (lines_index,) = nlte_data.lines_idx[excitation_species]

        try:
            j_blues_filtered = j_blues.iloc[lines_index]
        except AttributeError:
            j_blues_filtered = j_blues
        A_uls = nlte_data.A_uls[excitation_species]
        B_uls = nlte_data.B_uls[excitation_species]
        B_lus = nlte_data.B_lus[excitation_species]
        r_lu_index = lnu * number_of_levels + lnl
        r_ul_index = lnl * number_of_levels + lnu
        r_ul_matrix = np.zeros(
            (number_of_levels, number_of_levels, number_of_shells),
            dtype=np.float64,
        )
        r_ul_matrix_reshaped = r_ul_matrix.reshape(
            (number_of_levels**2, number_of_shells)
        )
        r_ul_matrix_reshaped[r_ul_index] = (
            A_uls[np.newaxis].T + B_uls[np.newaxis].T * j_blues_filtered
        )
        r_lu_matrix = np.zeros_like(r_ul_matrix)
        r_lu_matrix_reshaped = r_lu_matrix.reshape(
            (number_of_levels**2, number_of_shells)
        )
        r_lu_matrix_reshaped[r_lu_index] = (
            B_lus[np.newaxis].T * j_blues_filtered
        )
        return (
            lines_index,
            r_ul_index,
            r_ul_matrix,
            r_lu_index,
            r_lu_matrix,
        )
        # TODO: beta sobolev needs to be recalculated for each iteration, because it depends on number density

    def calculate_rate_matrix_nlte_exc(
        self,
        number_of_levels,
        nlte_data,
        j_blues,
        nlte_excitation_species,
    ):
        1 / 0
        return 0
