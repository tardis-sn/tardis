"""Calculate absorption probabilities using absorbing Markov chain theory.

This module implements the mathematical framework of absorbing Markov chains
to compute probabilities of photon absorption in each cell and the expected
number of steps before absorption from each source state.

References
----------
    Absorbing Markov chain theory: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import scipy


def create_absorbing_probs(
    transition_probabilities: pd.DataFrame,
    metadata: pd.DataFrame,
    max_workers: int = 1,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Calculate absorbing Markov chain probabilities and deactivation data.

    Computes the absorption probability matrices for each cell and extracts
    deactivation probabilities by solving the absorbing Markov chain system.
    The fundamental matrix is computed to determine absorption probabilities
    and expected number of steps to absorption from each state.

    WARNING: This is currently slower than it needs to be for many element sims.
    You should be able to slice out and process a subset of the probabilities
    based on some source state filtering. We have a draft PR open for this
    but it has not been finished.

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
    tuple[np.ndarray, pd.DataFrame]
        - absorbing_probability_matrix : np.ndarray
            Array of shape (num_cells, num_states, num_states) containing
            absorption probabilities from each state for each cell.
        - deactivating_probs : pd.DataFrame
            DataFrame of deactivation transition probabilities from
            absorption states.
    """
    num_cells = transition_probabilities.shape[1]

    internal_mask = metadata.transition_type >= 0
    internal_jump_probs = transition_probabilities[internal_mask]
    num_states = int(metadata.source_level_idx.max()) + 1

    absorbing_probability_matrix = np.zeros((num_cells, num_states, num_states))
    # Josh: The expected steps calculation is another linear algebra solve. We don't need
    # it for the MC calculation, so I'm removing it, but leaving it commented in case it is
    # useful for diagnostic purposes in the future.
    # expected_steps_in_cells_from_states = np.zeros((num_cells, num_states))

    rows = metadata[internal_mask].source_level_idx.values
    cols = metadata[internal_mask].destination_level_idx.values
    identity_matrix = scipy.sparse.identity(num_states, format="csc")

    def solve_cell(cell: int) -> tuple[np.ndarray, pd.DataFrame]:
        """Solve the absorbing Markov chain system for one cell.

        Parameters
        ----------
        cell : int
            Cell index in the transition probability table.

        Returns
        -------
        tuple[np.ndarray, pd.DataFrame]
            Full state-by-state absorbing probability matrix, and deactivating 
            probability DataFrame.
        """
        # In each cell, solve for absorbing markov chain probability
        # Follows math https://en.wikipedia.org/wiki/Absorbing_Markov_chain
        vals = internal_jump_probs.iloc[:, cell].to_numpy()

        internal_jump_matrix = scipy.sparse.coo_matrix(
            (vals, (rows, cols)), shape=(num_states, num_states)
        )

        identity_minus_Q = identity_matrix - internal_jump_matrix
        csc_N = identity_minus_Q.tocsc()

        # expected_steps = np.asarray(
        #     scipy.sparse.linalg.spsolve(csc_N, np.ones(num_states))
        # ).flatten()
        # expected_steps_in_cells_from_states[cell] = expected_steps

        deactivation_row = np.asarray(
            internal_jump_matrix.sum(axis=1)
        ).flatten()
        rhs_diagonal = 1 - deactivation_row
        active_rhs_columns = np.flatnonzero(rhs_diagonal)
        if len(active_rhs_columns) == 0:
            return cell, active_rhs_columns, np.empty((num_states, 0))

        # Solve (I - Q) * X = diag(1 - deactivation_row) only for columns
        # whose right-hand side is nonzero, then scatter back to the full
        # state-by-state matrix expected by transport.
        rhs = np.zeros((num_states, len(active_rhs_columns)))
        rhs[active_rhs_columns, np.arange(len(active_rhs_columns))] = (
            rhs_diagonal[active_rhs_columns]
        )
        solved_columns = np.asarray(
            scipy.sparse.linalg.splu(csc_N, permc_spec="MMD_AT_PLUS_A").solve(
                rhs
            )
        )
        if solved_columns.ndim == 1:
            solved_columns = solved_columns[:, np.newaxis]
        return cell, active_rhs_columns, solved_columns

    if max_workers > 1 and num_cells > 1:
        with ThreadPoolExecutor(
            max_workers=min(max_workers, num_cells)
        ) as executor:
            cell_solutions = executor.map(solve_cell, range(num_cells))
            for cell, active_rhs_columns, solved_columns in cell_solutions:
                absorbing_probability_matrix[cell][:, active_rhs_columns] = (
                    solved_columns
                )
    else:
        for cell in range(num_cells):
            _, active_rhs_columns, solved_columns = solve_cell(cell)
            absorbing_probability_matrix[cell][:, active_rhs_columns] = (
                solved_columns
            )

    deactivating_probs = transition_probabilities.copy()
    deactivating_probs[~internal_mask] = transition_probabilities[
        ~internal_mask
    ]

    return (
        absorbing_probability_matrix,
        deactivating_probs,
    )
