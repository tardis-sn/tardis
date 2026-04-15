import numpy as np
import pandas as pd


def create_absorbing_probs(
    transition_probabilities: pd.DataFrame, metadata: pd.DataFrame
):
    num_cells = transition_probabilities.shape[1]
    probs = transition_probabilities.copy()

    source_dest_pairs = list(
        zip(
            metadata.source_level_idx.values,
            metadata.destination_level_idx.values,
        )
    )
    probs["source_dest"] = source_dest_pairs

    num_states = len(metadata.source.unique())
    source_dest_index = pd.MultiIndex.from_product(
        [range(num_states), range(num_states)]
    )
    internal_mask = metadata.transition_type >= 0
    internal_jump_probs = probs[internal_mask]

    absorbing_probability_matrix = np.zeros((num_cells, num_states, num_states))
    expected_steps_in_cells_from_states = np.zeros((num_cells, num_states))

    for cell in range(num_cells):
        # In each cell, solve for absorbing markov chain probability
        # Follows math https://en.wikipedia.org/wiki/Absorbing_Markov_chain
        internal_jump_matrix = (
            internal_jump_probs.groupby("source_dest")
            .sum()
            .reindex(source_dest_index)
            .loc[source_dest_index, cell]
            .fillna(0)
            .values.reshape((num_states, num_states))
        )
        inv_matrix_fundamental_N = (
            np.identity(internal_jump_matrix.shape[0]) - internal_jump_matrix
        )
        fundamental_matrix_N = np.linalg.inv(inv_matrix_fundamental_N)

        expected_steps = fundamental_matrix_N.sum(axis=1)
        expected_steps_in_cells_from_states[cell] = expected_steps

        prob_deactivation_matrix = np.diag(1 - internal_jump_matrix.sum(axis=1))
        absorbing_probability_matrix[cell] = np.dot(
            fundamental_matrix_N, prob_deactivation_matrix
        )
    deactivating_probs = probs[~internal_mask].drop(columns="source_dest")
    deactivating_metadata = metadata[~internal_mask]

    return (
        absorbing_probability_matrix,
        deactivating_probs,
        deactivating_metadata,
    )
