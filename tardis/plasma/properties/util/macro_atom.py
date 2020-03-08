import numpy as np
from astropy import constants
from numba import njit, prange

h_cgs = constants.h.cgs.value
c = constants.c.cgs.value
kb = constants.k_B.cgs.value
inv_c2 = 1 / (c ** 2)

njit_dict = {'fastmath': False, 'parallel': False}

@njit(**njit_dict)
def calculate_transition_probabilities(
        transition_probability_coef,
        beta_sobolev, j_blues,
        stimulated_emission_factor,
        transition_type,
        lines_idx,
        block_references,
        transition_probabilities):

    norm_factor = np.zeros(transition_probabilities.shape[1])

    for i in prange(transition_probabilities.shape[0]):
        line_idx = lines_idx[i]
        for j in prange(transition_probabilities.shape[1]):
            transition_probabilities[i,j] = (transition_probability_coef[i] * beta_sobolev[line_idx,j])
        if transition_type[i] == 1:
            for j in prange(transition_probabilities.shape[1]):
                transition_probabilities[i,j] *= stimulated_emission_factor[line_idx,j] * j_blues[line_idx,j]

    for i in prange(block_references.shape[0] - 1):
        for k in prange(transition_probabilities.shape[1]):
            norm_factor[k] = 0.0
        for j in prange(block_references[i], block_references[i + 1]):
            for k in prange(transition_probabilities.shape[1]):
                norm_factor[k] += transition_probabilities[j,k]
        for k in prange(transition_probabilities.shape[1]):
            if norm_factor[k] != 0.0:
                norm_factor[k] = 1 / norm_factor[k]
            else:
                norm_factor[k] = 1.0
        for j in prange(block_references[i], block_references[i + 1]):
            for k in prange(0, transition_probabilities.shape[1]):
                transition_probabilities[j,k] *= norm_factor[k]


