import pandas as pd
import numpy as np
import pytest
from numpy.testing import assert_allclose


from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTERateEquationSolver,
)
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData


def test_prepare_bound_bound_rate_matrix(
    nlte_atomic_dataset,
):
    """
    Using a simple case of nlte_exc for HI, checks if prepare_bound_bound_rate_matrix generates the correct data.
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
    simple_j_blues_matrix = np.linspace(
        0.1, 1.0, copy_atomic_dataset.lines.index.size
    )
    simple_j_blues = pd.DataFrame(
        simple_j_blues_matrix,
        index=copy_atomic_dataset.lines.index,
        columns=["0"],
    )
    simple_number_of_levels = copy_atomic_dataset.levels.energy.loc[
        simple_excitation_species[0]
    ].count()
    (
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
    ) = NLTERateEquationSolver.prepare_r_uls_r_lus(
        simple_number_of_levels,
        simple_number_of_shells,
        simple_j_blues,
        simple_excitation_species[0],
        simple_nlte_data,
    )
    simple_beta_sobolev_matrix = np.linspace(
        2.5, 4.7, copy_atomic_dataset.lines.index.size
    )
    simple_beta_sobolev = pd.DataFrame(
        simple_beta_sobolev_matrix,
        index=copy_atomic_dataset.lines.index,
        columns=["0"],
    )
    actual_rate_matrix = NLTERateEquationSolver.prepare_bound_bound_rate_matrix(
        simple_number_of_levels,
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
        simple_beta_sobolev,
    )
    desired_rate_matrix = [
        [
            [-4.41229578e10],
            [1.09803977e10],
            [8.87031593e08],
            [1.83520728e08],
            [5.71742068e07],
        ],
        [
            [3.54409576e10],
            [-3.64473689e11],
            [1.32571818e11],
            [1.04228424e10],
            [2.27047121e09],
        ],
        [
            [5.71505717e09],
            [2.97836962e11],
            [-1.41199954e12],
            [6.39719360e11],
            [5.08905487e10],
        ],
        [
            [2.00818482e09],
            [4.15382182e10],
            [1.13720607e12],
            [-3.74593514e12],
            [1.98120242e12],
        ],
        [
            [9.58758249e08],
            [1.41181112e10],
            [1.41334621e11],
            [3.09560941e12],
            [-2.03442061e12],
        ],
    ]
    # if this test fails the first thing to check is if the reshape in the
    # methods made a view or a copy. If it's a copy rewrite the function.
    assert_allclose(
        desired_rate_matrix,
        np.array(actual_rate_matrix),
        rtol=1e-6,
    )


@pytest.mark.parametrize(
    [
        "coll_exc_coeff_values",
        "coll_deexc_coeff_values",
        "number_of_levels",
        "desired_coeff_matrix",
    ],
    [
        (
            [1, -2, 3],
            [4, 9, 10],
            3,
            [[1.0, 4.0, 9.0], [1.0, -7.0, 10.0], [-2.0, 3.0, -19.0]],
        ),
        (
            [0.21, 0.045, 0.1234],
            [0.7865, 0.987, 0.00123],
            3,
            [
                [-0.255, 0.7865, 0.987],
                [0.21, -0.9099, 0.00123],
                [0.045, 0.1234, -0.98823],
            ],
        ),
    ],
)
def test_coll_exc_deexc_matrix(
    coll_exc_coeff_values,
    coll_deexc_coeff_values,
    number_of_levels,
    desired_coeff_matrix,
):
    """
    Checks the NLTERateEquationSolver.create_coll_exc_deexc_matrix for simple values of species with 3 levels.
    NOTE: Values used for testing are not physical.
    """
    index = pd.MultiIndex.from_tuples(
        [(0, 1), (0, 2), (1, 2)],
        names=["level_number_lower", "level_number_upper"],
    )
    exc_coeff = pd.DataFrame(coll_exc_coeff_values, index=index)
    deexc_coeff = pd.DataFrame(coll_deexc_coeff_values, index=index)
    obtained_coeff_matrix = NLTERateEquationSolver.create_coll_exc_deexc_matrix(
        exc_coeff, deexc_coeff, number_of_levels
    )
    assert_allclose(obtained_coeff_matrix, desired_coeff_matrix)
