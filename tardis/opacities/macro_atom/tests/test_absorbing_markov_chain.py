import numpy as np
import pandas as pd
import scipy

from tardis.opacities.macro_atom.absorbing_markov_chain import (
    create_absorbing_probs,
)


def _create_absorbing_probs_reference(
    transition_probabilities: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[np.ndarray, pd.DataFrame]:
    num_cells = transition_probabilities.shape[1]
    num_states = len(metadata.source.unique())
    internal_mask = metadata.transition_type >= 0
    internal_jump_probs = transition_probabilities[internal_mask]

    absorbing_probability_matrix = np.zeros((num_cells, num_states, num_states))
    rows = metadata[internal_mask].source_level_idx.values
    cols = metadata[internal_mask].destination_level_idx.values
    for cell in range(num_cells):
        internal_jump_matrix = scipy.sparse.coo_matrix(
            (internal_jump_probs.iloc[:, cell].to_numpy(), (rows, cols)),
            shape=(num_states, num_states),
        )
        identity_minus_q = (
            scipy.sparse.identity(num_states) - internal_jump_matrix
        ).tocsc()
        deactivation_row = np.asarray(
            internal_jump_matrix.sum(axis=1)
        ).flatten()
        rhs = np.diag(1 - deactivation_row)
        absorbing_probability_matrix[cell] = scipy.sparse.linalg.spsolve(
            identity_minus_q, rhs
        )

    deactivating_probs = transition_probabilities.copy()
    deactivating_probs[internal_mask] *= 0
    return absorbing_probability_matrix, deactivating_probs


def test_create_absorbing_probs_matches_full_diagonal_solve() -> None:
    metadata = pd.DataFrame(
        {
            "source": [
                (1, 1, 0),
                (1, 1, 0),
                (1, 1, 1),
                (1, 1, 1),
                ("k", -99, -99),
                ("k", -99, -99),
            ],
            "destination": [
                (1, 1, 1),
                ("k", -99, -99),
                (1, 1, 0),
                ("k", -99, -99),
                (1, 1, 0),
                ("k", -99, -99),
            ],
            "transition_type": [1, -1, 1, -1, 1, -1],
            "source_level_idx": [0, 0, 1, 1, 2, 2],
            "destination_level_idx": [1, 2, 0, 2, 0, 2],
        }
    )
    transition_probabilities = pd.DataFrame(
        np.array(
            [
                [0.20, 0.15, 0.10],
                [0.80, 0.85, 0.90],
                [0.35, 0.25, 0.45],
                [0.65, 0.75, 0.55],
                [0.00, 0.00, 0.00],
                [1.00, 1.00, 1.00],
            ]
        )
    )

    actual_matrix, actual_deactivating_probs = create_absorbing_probs(
        transition_probabilities, metadata
    )
    expected_matrix, expected_deactivating_probs = (
        _create_absorbing_probs_reference(transition_probabilities, metadata)
    )

    np.testing.assert_allclose(actual_matrix, expected_matrix, rtol=1e-14)
    pd.testing.assert_frame_equal(
        actual_deactivating_probs, expected_deactivating_probs
    )
