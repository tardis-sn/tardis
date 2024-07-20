import pandas as pd
import numpy as np

from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTERateEquationsSolver,
)


def create_nlte_excitation_rate_matrix_block(
    atomic_number,
    ion_number,
    gamma,
    alpha_sp,
    alpha_stim,
    coll_ion_coeff,
    coll_recomb_coeff,
    coll_exc_coeff,
    coll_deexc_coeff,
    electron_density,
    nlte_excitation_species,
    j_blues,
    nlte_data,
    beta_sobolev,
):
    number_of_levels = gamma.loc[(atomic_number, ion_number)].size
    number_of_shells = 1
    ionization_coeff = gamma + coll_ion_coeff * electron_density
    recombination_coeff = (
        alpha_sp + alpha_stim
    ) * electron_density + coll_recomb_coeff * electron_density**2
    diagonal_block = np.diag(
        -ionization_coeff.loc[(atomic_number, ion_number)]
        - recombination_coeff.loc[(atomic_number, ion_number)]
    )
    ion_recomb_block = np.zeros(
        (diagonal_block.shape[0] + 1, diagonal_block.shape[1] + 1)
    )
    if ion_number > 0:
        try:
            ion_recomb_block[1:-1, 0] = ionization_coeff.loc[
                (atomic_number, ion_number - 1)
            ]
        except:
            ion_recomb_block[1:-1, 0] = 0
        try:
            ionization_coeff[0, 1:-1] = recombination_coeff.loc[
                (atomic_number, ion_number - 1)
            ]
        except:
            ionization_coeff[0, 1:-1] = 0
    if ion_number < atomic_number:
        try:
            ion_recomb_block[-1, 1:-1] = ionization_coeff.loc[
                (atomic_number, ion_number)
            ]
        except:
            ion_recomb_block[-1, 1:-1] = 0
        try:
            ion_recomb_block[1:-1, -1] = recombination_coeff
        except:
            ion_recomb_block[1:-1, -1] = 0

    coll_exc_deexc_block = (
        NLTERateEquationsSolver.create_coll_exc_deexc_block(
            coll_exc_coeff, coll_deexc_coeff, number_of_levels
        )
        * electron_density
    )

    (
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
    ) = NLTERateEquationsSolver.prepare_r_uls_r_lus(
        number_of_levels,
        number_of_shells,
        j_blues,
        nlte_excitation_species,
        nlte_data,
    )
    bound_bound_block = NLTERateEquationsSolver.prepare_bound_bound_rate_matrix(
        number_of_levels,
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
        beta_sobolev,
    )
    excitation_block = bound_bound_block + coll_exc_deexc_block
    rate_matrix_block = ion_recomb_block
    rate_matrix_block[1:-1, 1:-1] = excitation_block
    return rate_matrix_block
