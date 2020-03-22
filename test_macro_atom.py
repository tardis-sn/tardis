from tardis.plasma.properties.util import macro_atom, macro_atom_with_numba, macro_atom_with_numpy
# Check for all three scripts one-by-one
import numpy as np
import time

n = 5000
transition_probability_coef = np.random.rand(n)
beta_sobolev = np.random.rand(n, n)
j_blues = np.random.rand(n, n)
stimulated_emission_factor = np.random.rand(n, n)
transition_type = np.int64(np.random.rand(n))
lines_idx = np.int64(np.random.rand(n))
block_references = np.int64(np.random.rand(n))
transition_probabilities = np.random.rand(n, n)

start = time.time()
# Replace with different script
macro_atom_with_numpy.calculate_transition_probabilities(
        transition_probability_coef,
        beta_sobolev, j_blues,
        stimulated_emission_factor,
        transition_type,
        lines_idx,
        block_references,
        transition_probabilities)
t1 = time.time()-start

start = time.time()
# Replace with different script
macro_atom_with_numpy.calculate_transition_probabilities(
        transition_probability_coef,
        beta_sobolev, j_blues,
        stimulated_emission_factor,
        transition_type,
        lines_idx,
        block_references,
        transition_probabilities)
t2 = time.time()-start

print("Time taken (with compilation): {}".format(t1))
print("Time taken (after compilation): {}".format(t2))

