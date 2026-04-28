"""Calculate absorption probabilities using absorbing Markov chain theory.

This module implements the mathematical framework of absorbing Markov chains
to compute probabilities of photon absorption in each cell and the expected
number of steps before absorption from each source state.

References:
    Absorbing Markov chain theory: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
"""

import numpy as np
import pandas as pd


def create_absorbing_probs(
    transition_probabilities: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Calculate absorbing Markov chain probabilities and deactivation data.

    Computes the absorption probability matrices for each cell and extracts
    deactivation probabilities by solving the absorbing Markov chain system.
    The fundamental matrix is computed to determine absorption probabilities
    and expected number of steps to absorption from each state.

    Parameters
    ----------
    transition_probabilities : pd.DataFrame
        DataFrame with transition probabilities between states for each cell.
        Rows represent transitions, columns represent cells.
    metadata : pd.DataFrame
        Metadata about transitions including source_level_idx,
        destination_level_idx, and transition_type. Values >= 0 for
        transition_type indicate internal transitions transitions that are not accompanied
        by reemision of the r-packet, and an end to the interaction handler chain.
        Negative values indicate deactivation, or reemitting from an absorbing state.

    Returns
    -------
    tuple[np.ndarray, pd.DataFrame, pd.DataFrame]
        - absorbing_probability_matrix : np.ndarray
            Array of shape (num_cells, num_states, num_states) containing
            absorption probabilities from each state for each cell.
        - deactivating_probs : pd.DataFrame
            DataFrame of deactivation transition probabilities from
            absorption states.
        - deactivating_metadata : pd.DataFrame
            Metadata for deactivation transitions.
    """
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
