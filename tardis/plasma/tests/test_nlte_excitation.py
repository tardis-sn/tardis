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
    simple_j_blues = pd.DataFrame(
        0.5, index=copy_atomic_dataset.lines.index, columns=["0"]
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
    simple_beta_sobolev = pd.DataFrame(
        0.8, index=copy_atomic_dataset.lines.index, columns=["0"]
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
            [-4.22105726e09],
            [1.22518171e09],
            [1.05118978e08],
            [2.20562167e07],
            [6.85522197e06],
        ],
        [
            [3.39792096e09],
            [-3.46604468e10],
            [1.25856407e10],
            [9.80445669e08],
            [2.10361778e08],
        ],
        [
            [5.44879921e08],
            [2.82383118e10],
            [-1.31698271e11],
            [5.96565953e10],
            [4.66884472e09],
        ],
        [
            [1.89343842e08],
            [3.89484475e09],
            [1.06043389e11],
            [-3.43063731e11],
            [1.80741125e11],
        ],
        [
            [8.89125366e07],
            [1.30210860e09],
            [1.29641220e10],
            [2.82404633e11],
            [-1.85627187e11],
        ],
    ]
    breakpoint()

    assert_allclose(
        desired_rate_matrix,
        np.array(actual_rate_matrix),
        rtol=1e-6,
    )
