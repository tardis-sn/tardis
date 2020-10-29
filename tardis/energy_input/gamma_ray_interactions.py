import numpy as np

from tardis.energy_input.util import kappa_calculation, klein_nishina, euler_rodrigues, compton_theta_distribution
from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu

def compton_scatter(gamma_ray):
    theta_angles, theta_distribution = compton_theta_distribution(gamma_ray.energy)
    
    z = np.random.random()
    
    #sample new random theta direction
    random_vector = get_random_mu()
    
    #REWRITE FOR 1D
    perpendicular_vector_x = gamma_ray.direction.y * random_vector.z - gamma_ray.direction.z * random_vector.y
    perpendicular_vector_y = gamma_ray.direction.z * random_vector.x - gamma_ray.direction.x * random_vector.z
    perpendicular_vector_z = gamma_ray.direction.x * random_vector.y - gamma_ray.direction.y * random_vector.x
    
    perpendicular_vector = np.array([perpendicular_vector_x, perpendicular_vector_y, perpendicular_vector_z])
    
    #get Compton scattering angle
    compton_angle = theta_angles[np.searchsorted(theta_distribution, z)]
    
    #rotate to match
    rotation_matrix = euler_rodrigues(compton_angle, gamma_ray.direction)
    resulting_direction = np.dot(rotation_matrix, perpendicular_vector)
    
    gamma_ray.direction.mu = resulting_direction[2]
    
    #Energy calculations
    new_energy = gamma_ray.energy / (1. + kappa_calculation(gamma_ray.energy) * (1. - gamma_ray.direction.mu))  
    lost_energy = gamma_ray.energy - new_energy   
    gamma_ray.energy = new_energy
    
    return lost_energy

def pair_creation(gamma_ray):
    direction_mu = get_random_mu()
    
    gamma_ray.energy = 511.0e3
    gamma_ray.direction = direction_mu

def photoabsorption(gamma_ray):
    gamma_ray.status = 'Absorbed'
    return gamma_ray.energy