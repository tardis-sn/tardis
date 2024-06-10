import numpy as np
from numba import njit

from tardis import constants as const
from tardis.transport.montecarlo import njit_dict


@njit(**njit_dict)
def fast_calculate_transition_probabilities(
    transition_probability_coef,
    beta_sobolev,
    j_blues,
    stimulated_emission_factor,
    transition_type,
    lines_idx,
    block_references,
    transition_probabilities,
    normalize,
):
    """
    Calculates transition probabilities for macro_atom interactions

    transition_probability_coef must be a 1D array
    transition_type, lines_idx, and block_references must be int-type arrays
    beta_sobolev, j_blues,stimulated_emission_factor, and transition_probabilities must be 2D array
    """
    norm_factor = np.zeros(transition_probabilities.shape[1])

    for i in range(transition_probabilities.shape[0]):
        line_idx = lines_idx[i]
        for j in range(transition_probabilities.shape[1]):
            transition_probabilities[i, j] = (
                transition_probability_coef[i] * beta_sobolev[line_idx, j]
            )
        if transition_type[i] == 1:
            for j in range(transition_probabilities.shape[1]):
                transition_probabilities[i, j] *= (
                    stimulated_emission_factor[line_idx, j]
                    * j_blues[line_idx, j]
                )

    if normalize:
        for i in range(block_references.shape[0] - 1):
            for k in range(transition_probabilities.shape[1]):
                norm_factor[k] = 0.0
            for j in range(block_references[i], block_references[i + 1]):
                for k in range(transition_probabilities.shape[1]):
                    norm_factor[k] += transition_probabilities[j, k]
            for k in range(transition_probabilities.shape[1]):
                if norm_factor[k] != 0.0:
                    norm_factor[k] = 1 / norm_factor[k]
                else:
                    norm_factor[k] = 1.0
            for j in range(block_references[i], block_references[i + 1]):
                for k in range(transition_probabilities.shape[1]):
                    transition_probabilities[j, k] *= norm_factor[k]
