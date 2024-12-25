import numpy as np
from numba import float64
from numba.experimental import jitclass

def initialize_j_blue_estimator(tau_sobolev_shape):

    j_blue_estimator = np.zeros(tau_sobolev_shape)
    return JBlueEstimator(j_blue_estimator)

@jitclass([
    ("j_blue_estimator", float64[:, :]),
])
class JBlueEstimator:
    def __init__(self, j_blue_estimator):
        self.j_blue_estimator = j_blue_estimator

    def increment(self, other):

        self.j_blue_estimator += other.j_blue_estimator 

if __name__ == "__main__":
    # Create a small example with 3x3 shape
    test_estimator = initialize_j_blue_estimator((3, 3))
    print("Initial j_blue_estimator:")
    print(test_estimator.j_blue_estimator)

    # Create another estimator and modify it to show increment
    other_estimator = initialize_j_blue_estimator((3, 3))
    other_estimator.j_blue_estimator[0, 0] = 1.0
    other_estimator.j_blue_estimator[1, 1] = 2.0

    # Increment the first estimator with the second one
    test_estimator.increment(other_estimator)
    print("\nAfter incrementing with another estimator:")
    print(test_estimator.j_blue_estimator)

    # Save the estimator to a .npy file
    output_file = "j_blue_estimator.npy"
    np.save(output_file, test_estimator.j_blue_estimator)
    print(f"\nSaved estimator to {output_file}") 