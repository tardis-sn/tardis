import pandas as pd
import numpy as np
from numpy.testing import assert_almost_equal


from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTERateEquationSolver,
)
from tardis.plasma.properties.nlte_excitation_data import NLTEExcitationData


def test_main_nlte_calculation_bound_bound(
    nlte_atomic_dataset,
):
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """
    simple_excitation_species = (1, 0)
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
    simple_t_electrons = [0.5646738847]
    simple_j_blues = pd.DataFrame(
        0.5, index=copy_atomic_dataset.lines.index, columns=["0"]
    )

    simple_beta_sobolev = pd.DataFrame(
        0.8, index=copy_atomic_dataset.lines.index, columns=["0"]
    )

    actual_rate_matrix = (
        NLTERateEquationSolver.main_nlte_calculation_bound_bound(
            copy_atomic_dataset.levels,
            simple_t_electrons,
            simple_j_blues,
            simple_beta_sobolev,
            simple_excitation_species,
            simple_nlte_data,
        )
    )
    desired_rate_matrix = [
        [
            [1.00000000e00],
            [1.00000000e00],
            [1.00000000e00],
            [1.00000000e00],
            [1.00000000e00],
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

    assert_almost_equal(
        desired_rate_matrix / np.array(actual_rate_matrix),
        np.ones_like(desired_rate_matrix),
        decimal=6,
    )
