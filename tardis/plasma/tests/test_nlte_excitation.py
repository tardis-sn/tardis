import pandas as pd
import numpy as np
from numpy.testing import assert_allclose


from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTERateEquationSolver,
)
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData


def test_prepare_bound_bound_rate_matrix(
    nlte_atomic_dataset,
):
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """
    simple_excitation_species = [(1, 0)]
    copy_atomic_dataset = nlte_atomic_dataset
    copy_atomic_dataset.levels = nlte_atomic_dataset.levels[
        nlte_atomic_dataset.levels.index.get_level_values("level_number") < 5
    ]
    lines_filtered = nlte_atomic_dataset.lines[
        nlte_atomic_dataset.lines.index.get_level_values("level_number_lower")
        < 5
    ]
    copy_atomic_dataset.lines = lines_filtered[
        lines_filtered.index.get_level_values("level_number_upper") < 5
    ]
    simple_nlte_data = NLTEExcitationData(
        copy_atomic_dataset.lines, simple_excitation_species
    )
    simple_number_of_shells = 1
    simple_j_blues_matrix = np.linspace(0.1, 1., copy_atomic_dataset.lines.index.size)
    simple_j_blues = pd.DataFrame(
        simple_j_blues_matrix, index=copy_atomic_dataset.lines.index, columns=["0"]
    )

    (
        lines_index,
        number_of_levels,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
    ) = NLTERateEquationSolver.prepare_r_uls_rlus(
        copy_atomic_dataset.levels,
        simple_number_of_shells,
        simple_j_blues,
        simple_excitation_species[0],
        simple_nlte_data,
    )
    simple_beta_sobolev_matrix = np.linspace(2.5, 4.7, copy_atomic_dataset.lines.index.size)
    simple_beta_sobolev = pd.DataFrame(
        simple_beta_sobolev_matrix, index=copy_atomic_dataset.lines.index, columns=["0"]
    )
    actual_rate_matrix = NLTERateEquationSolver.prepare_bound_bound_rate_matrix(
        number_of_levels,
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
        simple_beta_sobolev,
    )
    desired_rate_matrix = [
        [
            [-4.41229578e+10],
            [1.09803977e+10],
            [8.87031593e+08],
            [1.83520728e+08],
            [5.71742068e+07],
        ],
        [
            [3.54409576e+10],
            [-3.64473689e+11],
            [1.32571818e+11],
            [1.04228424e+10],
            [2.27047121e+09],
        ],
        [
            [5.71505717e+09],
            [2.97836962e+11],
            [-1.41199954e+12],
            [6.39719360e+11],
            [5.08905487e+10],
        ],
        [
            [2.00818482e+09],
            [4.15382182e+10],
            [1.13720607e+12],
            [-3.74593514e+12],
            [1.98120242e+12],
        ],
        [
            [9.58758249e+08],
            [1.41181112e+10],
            [1.41334621e+11],
            [3.09560941e+12],
            [-2.03442061e+12],
        ],
    ]
    assert_allclose(
        desired_rate_matrix,
        np.array(actual_rate_matrix),
        rtol=1e-6,
    )
