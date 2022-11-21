import pytest
import numpy as np
import pandas as pd
from numpy.testing import assert_almost_equal
from tardis.plasma.properties import NLTERateEquationSolver


def test_rate_matrix():
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """

    simple_index_nlte_ion = pd.MultiIndex.from_tuples(
        [(1, 0), (2, 1)], names=("atomic_number", "ion_number")
    )
    simple_index_lte_ion = pd.MultiIndex.from_tuples(
        [(1, 1), (2, 1), (2, 2)], names=("atomic_number", "ion_number")
    )
    simple_rate_matrix_index = pd.MultiIndex.from_tuples(
        [
            (1, 0, "nlte_ion"),
            (1, 1, "lte_ion"),
            (2, 0, "lte_ion"),
            (2, 1, "nlte_ion"),
            (2, 2, "lte_ion"),
            ("n_e", "n_e", "n_e"),
        ],
        names=("atomic_number", "ion_number", "level_number"),
    )

    simple_photo_ion_rates = [0.03464792, 0.68099508]
    simple_rad_recomb_coeff = [0.43303813, 0.66140309]
    simple_col_ion_coeff = [0.19351674, 0.69214007]
    simple_col_recomb_coeff = [0.06402515, 0.29785023]
    simple_phi = [0.18936306, 0.15726292, 0.79851244]

    simple_photo_ion_rates = pd.DataFrame(
        simple_photo_ion_rates, index=simple_index_nlte_ion
    )
    simple_rad_recomb_coeff = pd.DataFrame(
        simple_rad_recomb_coeff, index=simple_index_nlte_ion
    )
    simple_col_ion_coeff = pd.DataFrame(
        simple_col_ion_coeff, index=simple_index_nlte_ion
    )
    simple_col_recomb_coeff = pd.DataFrame(
        simple_col_recomb_coeff, index=simple_index_nlte_ion
    )
    simple_phi = pd.DataFrame(simple_phi, index=simple_index_lte_ion)

    simple_electron_density = 0.2219604493076

    actual_rate_matrix = NLTERateEquationSolver.calculate_rate_matrix(
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_photo_ion_rates,
        simple_rad_recomb_coeff,
        simple_col_ion_coeff,
        simple_col_recomb_coeff,
    )
    desired_rate_matrix = [
        [-0.077601, 0.099272, 0.000000, 0.000000, 0.000000, 0.0],
        [1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.0],
        [0.000000, 0.000000, -0.157263, 0.221960, 0.000000, 0.0],
        [0.000000, 0.000000, 0.000000, -0.834623, 0.161479, 0.0],
        [0.000000, 0.000000, 1.000000, 1.000000, 1.000000, 0.0],
        [0.000000, 1.000000, 0.000000, 1.000000, 2.000000, -1.0],
    ]

    assert_almost_equal(
        desired_rate_matrix, np.array(actual_rate_matrix), decimal=6
    )


def test_jacobian_matrix():
    simple_index_nlte_ion = pd.MultiIndex.from_tuples(
        [(1, 0), (2, 1)], names=("atomic_number", "ion_number")
    )
    simple_index_lte_ion = pd.MultiIndex.from_tuples(
        [(1, 1), (2, 1), (2, 2)], names=("atomic_number", "ion_number")
    )
    simple_rate_matrix_index = pd.MultiIndex.from_tuples(
        [
            (1, 0, "nlte_ion"),
            (1, 1, "lte_ion"),
            (2, 0, "lte_ion"),
            (2, 1, "nlte_ion"),
            (2, 2, "lte_ion"),
            ("n_e", "n_e", "n_e"),
        ],
        names=("atomic_number", "ion_number", "level_number"),
    )

    simple_photo_ion_rates = [0.03464792, 0.68099508]
    simple_rad_recomb_coeff = [0.43303813, 0.66140309]
    simple_col_ion_coeff = [0.19351674, 0.69214007]
    simple_col_recomb_coeff = [0.06402515, 0.29785023]
    simple_phi = [0.18936306, 0.15726292, 0.79851244]

    simple_photo_ion_rates = pd.DataFrame(
        simple_photo_ion_rates, index=simple_index_nlte_ion
    )
    simple_rad_recomb_coeff = pd.DataFrame(
        simple_rad_recomb_coeff, index=simple_index_nlte_ion
    )
    simple_col_ion_coeff = pd.DataFrame(
        simple_col_ion_coeff, index=simple_index_nlte_ion
    )
    simple_col_recomb_coeff = pd.DataFrame(
        simple_col_recomb_coeff, index=simple_index_nlte_ion
    )
    simple_phi = pd.DataFrame(simple_phi, index=simple_index_lte_ion)

    simple_electron_density = 0.2219604493076

    initial_guess = [
        0.7192433675307516,
        0.8101666197902874,
        0.7171853313284426,
        0.040220760173800496,
        0.2878574499274399,
        simple_electron_density,
    ]
    simple_rate_matrix = NLTERateEquationSolver.calculate_rate_matrix(
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_photo_ion_rates,
        simple_rad_recomb_coeff,
        simple_col_ion_coeff,
        simple_col_recomb_coeff,
    )

    actual_jacobian_matrix = NLTERateEquationSolver.jacobian_matrix(
        initial_guess,
        simple_rate_matrix,
        simple_rate_matrix_index,
        simple_rad_recomb_coeff,
        simple_col_ion_coeff,
        simple_col_recomb_coeff,
    )

    desired_jacobian_matrix = [
        [-0.07760098, 0.09927163, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.15726292, 0.22196045, 0.0, 0.04022076],
        [0.0, 0.0, 0.0, -0.8346228, 0.16147935, 0.0],
        [0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 2.0, -1.0],
    ]
    assert_almost_equal(actual_jacobian_matrix, desired_jacobian_matrix)
