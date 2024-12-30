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
    # Create multiple test cases with different shapes
    small_estimator = initialize_j_blue_estimator((5, 10))
    large_estimator = initialize_j_blue_estimator((20, 40))
    
    # Simulate different packet interactions
    packet_energies = [1.0e-10, 2.0e-10, 5.0e-10]
    packet_nus = [2.0e15, 3.0e15, 4.0e15]
    
    # Add some interactions to small estimator
    for energy, nu in zip(packet_energies, packet_nus):
        small_estimator.j_blue_estimator[2, 5] += energy / nu
        small_estimator.j_blue_estimator[3, 7] += energy / (2 * nu)
    
    # Add some interactions to large estimator
    for energy, nu in zip(packet_energies, packet_nus):
        large_estimator.j_blue_estimator[10, 20] += energy / nu
        large_estimator.j_blue_estimator[15, 30] += energy / (3 * nu)
    
    # Print raw values
    print("Small estimator interactions:")
    print(f"Position (2,5): {small_estimator.j_blue_estimator[2, 5]}")
    print(f"Position (3,7): {small_estimator.j_blue_estimator[3, 7]}")
    
    print("\nLarge estimator interactions:")
    print(f"Position (10,20): {large_estimator.j_blue_estimator[10, 20]}")
    print(f"Position (15,30): {large_estimator.j_blue_estimator[15, 30]}")
    
    # Test estimator increment
    small_estimator.increment(small_estimator)
    print("\nAfter incrementing small estimator with itself:")
    print(f"Position (2,5): {small_estimator.j_blue_estimator[2, 5]}")
    
    # Simulate different normalizations
    time_explosions = [1.0e5, 2.0e5]
    time_simulations = [1.0e4, 2.0e4]
    volumes = [1.0e45, 2.0e45]
    
    print("\nTesting different normalizations:")
    for t_explosion, t_simulation, vol in zip(time_explosions, time_simulations, volumes):
        norm_factor = 3e10 * t_explosion / (4 * np.pi * t_simulation * vol)
        normalized = small_estimator.j_blue_estimator * norm_factor
        print(f"\nNormalization with:")
        print(f"t_explosion={t_explosion}, t_simulation={t_simulation}, volume={vol}")
        print(f"Position (2,5): {normalized[2, 5]}")
    
    # Save examples to files
    np.save("j_blue_small_raw.npy", small_estimator.j_blue_estimator)
    np.save("j_blue_large_raw.npy", large_estimator.j_blue_estimator) 