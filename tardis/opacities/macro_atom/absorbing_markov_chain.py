"""Calculate absorption probabilities using absorbing Markov chain theory.

This module implements the mathematical framework of absorbing Markov chains
to compute probabilities of photon absorption in each cell and the expected
number of steps before absorption from each source state.

References:
    Absorbing Markov chain theory: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
"""

import numpy as np
import pandas as pd
import scipy


def create_absorbing_probs(
    transition_probabilities: pd.DataFrame,
    metadata: pd.DataFrame,
    selected_species=None,
) -> tuple[np.ndarray, pd.DataFrame, dict]:
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
        by reemission of the r-packet, and an end to the interaction handler chain.
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
        - source_idx_to_absorbing_matrix_idx_map : dict
            Dictionary describing what each level of the absorbing probability matrix
            maps to.
    """
    num_cells = transition_probabilities.shape[1]

    internal_mask = metadata.transition_type >= 0

    selected_species = [(1, 0)]  # TODO: FIX TO PASS CORRECTLY FROM CONFIG
    if selected_species:
        selected_species = (
            list(selected_species)
            + [
                ("i", -99),
                # ("k", -99),
            ]
        )  # We also want the i and k blocks to be include in our selection, or we won't process through the thermal pool

        selected_species_mask = np.zeros_like(
            metadata.source.values, dtype=bool
        )

        sliced_source_metadata_col = metadata.source.apply(lambda x: x[:2])
        sliced_dest_metadata_col = metadata.destination.apply(lambda x: x[:2])
        # build a mask to grab all transitions that include a selected species and k and i block
        # JOSH: Will need to revisit with i-block handling, maybe

        # This creates the map we need to point to the correct abs markov spot

        for species in selected_species:
            species_source_mask = sliced_source_metadata_col == species
            selected_species_mask = np.logical_or(
                selected_species_mask, species_source_mask
            )
            species_dest_mask = sliced_dest_metadata_col == species
            selected_species_mask = np.logical_or(
                selected_species_mask, species_dest_mask
            )

        internal_selected_mask = np.logical_and(
            selected_species_mask, internal_mask
        )  # we want the internal jumps of the selected species

        source_idx_to_absorbing_matrix_idx_map = {
            np.int64(source_level_idx): np.int64(abs_block_idx)
            for abs_block_idx, source_level_idx in enumerate(
                metadata[internal_selected_mask].source_level_idx.unique()
            )
        }  # This will tell you which block of the matrix you should enter, and which exit you go to because the matrix is symmetric

    else:
        internal_selected_mask = internal_mask
        source_idx_to_absorbing_matrix_idx_map = {}

    num_states = len(metadata[internal_selected_mask].source.unique())

    internal_jump_probs = transition_probabilities[internal_selected_mask]

    absorbing_probability_matrix = np.zeros((num_cells, num_states, num_states))
    # Josh: The expected steps calculation is another linear algebra solve. We don't need
    # it for the MC calculation, so I'm removing it, but leaving it commented in case it is
    # useful for diagnostic purposes in the future.
    # expected_steps_in_cells_from_states = np.zeros((num_cells, num_states))

    rows = metadata[internal_selected_mask].source_level_idx.values
    cols = metadata[internal_selected_mask].destination_level_idx.values

    matrix_rows = [source_idx_to_absorbing_matrix_idx_map[x] for x in rows]
    matrix_cols = [source_idx_to_absorbing_matrix_idx_map[x] for x in cols]

    for cell in range(num_cells):
        print(f"starting cell {cell}")
        # In each cell, solve for absorbing markov chain probability
        # Follows math https://en.wikipedia.org/wiki/Absorbing_Markov_chain
        vals = internal_jump_probs[cell].values

        internal_jump_matrix = scipy.sparse.coo_matrix(
            (vals, (matrix_rows, matrix_cols)), shape=(num_states, num_states)
        )

        identity_minus_Q = (
            scipy.sparse.identity(internal_jump_matrix.shape[0])
            - internal_jump_matrix
        )
        csc_N = identity_minus_Q.tocsc()

        # expected_steps = np.asarray(
        #     scipy.sparse.linalg.spsolve(csc_N, np.ones(num_states))
        # ).flatten()
        # expected_steps_in_cells_from_states[cell] = expected_steps

        deactivation_row = np.asarray(
            internal_jump_matrix.sum(axis=1)
        ).flatten()
        # Solve (I - Q) * X = diag(1 - deactivation_row) directly
        # instead of computing inv(I - Q) explicitly
        rhs = np.diag(1 - deactivation_row)
        absorbing_probability_matrix[cell] = scipy.sparse.linalg.spsolve(
            csc_N, rhs
        )

    deactivating_probs = transition_probabilities.copy()
    deactivating_probs[internal_selected_mask] *= 0

    return (
        absorbing_probability_matrix,
        deactivating_probs,
        source_idx_to_absorbing_matrix_idx_map,
    )
