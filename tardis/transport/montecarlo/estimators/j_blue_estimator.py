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
    # Create example with realistic shape
    test_estimator = initialize_j_blue_estimator((10, 20))  # typical tau_sobolev shape
    
    # Simulate some updates
    packet_energy = 1.0e-10
    packet_nu = 2.0e15
    test_estimator.j_blue_estimator[5, 10] += packet_energy / packet_nu
    
    print("After simulated packet interaction:")
    print(test_estimator.j_blue_estimator[5, 10])
    
    # Simulate normalization
    time_explosion = 1.0e5
    time_simulation = 1.0e4
    volume = 1.0e45
    norm_factor = 3e10 * time_explosion / (4 * np.pi * time_simulation * volume)
    
    normalized = test_estimator.j_blue_estimator * norm_factor
    print("\nAfter normalization:")
    print(normalized[5, 10])
    
    np.save("j_blue_raw.npy", test_estimator.j_blue_estimator)
    np.save("j_blue_normalized.npy", normalized) 