import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData
from tardis.plasma.properties.nlte_rate_equation_solver import (
    prepare_r_uls_r_lus,
    prepare_bound_bound_rate_matrix,
    create_coll_exc_deexc_matrix,
)


def test_prepare_bound_bound_rate_matrix(nlte_atomic_dataset, regression_data):
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
    ) = prepare_r_uls_r_lus(
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
    actual_rate_matrix = prepare_bound_bound_rate_matrix(
        simple_number_of_levels,
        lines_index,
        r_ul_index,
        r_ul_matrix,
        r_lu_index,
        r_lu_matrix,
        simple_beta_sobolev,
    )
    # if this test fails the first thing to check is if the reshape in the
    # methods made a view or a copy. If it's a copy rewrite the function.
    # TODO: allow rtol=1e-6
    expected_rate_matrix = regression_data.sync_ndarray(actual_rate_matrix)
    npt.assert_allclose(actual_rate_matrix, expected_rate_matrix, rtol=1e-6)


@pytest.mark.parametrize(
    [
        "coll_exc_coeff_values",
        "coll_deexc_coeff_values",
        "number_of_levels",
    ],
    [
        (
            [1, -2, 3],
            [4, 9, 10],
            3,
        ),
        (
            [0.21, 0.045, 0.1234],
            [0.7865, 0.987, 0.00123],
            3,
        ),
    ],
)
def test_coll_exc_deexc_matrix(
    coll_exc_coeff_values,
    coll_deexc_coeff_values,
    number_of_levels,
    regression_data,
):
    """
    Checks the NLTERateEquationSolver.create_coll_exc_deexc_matrix for simple values of species with 3 levels.
    NOTE: Values used for testing are not physical.
    """
    index = pd.MultiIndex.from_tuples(
        [(0, 1), (0, 2), (1, 2)],
        names=["level_number_lower", "level_number_upper"],
    )
    exc_coeff = pd.Series(coll_exc_coeff_values, index=index)
    deexc_coeff = pd.Series(coll_deexc_coeff_values, index=index)
    obtained_coeff_matrix = create_coll_exc_deexc_matrix(
        exc_coeff, deexc_coeff, number_of_levels
    )
    expected_obtained_coeff_matrix = regression_data.sync_ndarray(
        obtained_coeff_matrix
    )
    npt.assert_allclose(expected_obtained_coeff_matrix, obtained_coeff_matrix)
