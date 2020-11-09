import numpy as np
from tardis.energy_input.util import quadratic

def calculate_distance_radial(gamma_ray, r_inner, r_outer):
    
    position_square = gamma_ray.location.x ** 2. + gamma_ray.location.y ** 2. + gamma_ray.location.z ** 2.
    position_direction = gamma_ray.direction.x * gamma_ray.location.x + \
                        gamma_ray.direction.y * gamma_ray.location.y + \
                        gamma_ray.direction.z * gamma_ray.location.z
    
    on_inner_wall = False
    on_outer_wall = False
    
    distances = []
    
    if np.abs(gamma_ray.location.r - r_inner) < 5:
        on_inner_wall = True
        
    if np.abs(gamma_ray.location.r - r_outer) < 5:
        on_outer_wall = True

    quadratic_b = 2. * position_direction
    quadratic_c = position_square
    
    quad_c_inner = quadratic_c - r_inner ** 2
    distance_1, distance_2 = quadratic(quadratic_b, quad_c_inner)
    
    if on_inner_wall:
        if np.abs(distance_1) < np.abs(distance_2):
            distances.append(distance_2)
        else:
            distances.append(distance_1)
    else:
        distances.append(distance_1)
        distances.append(distance_2)
    
    quad_c_outer = quadratic_c - r_outer ** 2
    distance_3, distance_4 = quadratic(quadratic_b, quad_c_outer)
    
    if on_outer_wall:
        if np.abs(distance_3) < np.abs(distance_4):
            distances.append(distance_4)
        else:
            distances.append(distance_3)
    else:
        distances.append(distance_3)
        distances.append(distance_4)
        
    distances = [item for item in distances if item >= 0]
    if len(distances) > 0:
        distance = min(distances, key=abs)
    else:
        distance = 0.0

    return distance

def distance_trace(gamma_ray, radii, total_opacity, distance_moved):
    
    if gamma_ray.shell < len(radii) - 1:
        distance_boundary = calculate_distance_radial(gamma_ray, radii[gamma_ray.shell], radii[packet.shell + 1])
    else:
        distance_boundary = 0.0
    
    z = np.random.random()
    tau = -np.log(z)
    
    distance_interaction = tau / total_opacity - distance_moved
    
    if distance_interaction < distance_boundary:
        return distance_interaction, distance_boundary, True
    else:
        return distance_interaction, distance_boundary, False

def move_gamma_ray(gamma_ray, distance):
    x_old = gamma_ray.location.x
    y_old = gamma_ray.location.y
    z_old = gamma_ray.location.z
    
    x_new = x_old + distance * gamma_ray.direction.x
    y_new = y_old + distance * gamma_ray.direction.y
    z_new = z_old + distance * gamma_ray.direction.z
    
    gamma_ray.location.r = np.sqrt(x_new ** 2. + y_new ** 2. + z_new ** 2.)
    gamma_ray.location.mu = z_new / gamma_ray.location.r
    
    return gamma_ray

def density_sampler(radii, mass_ratio):   
    z = np.random.random()   
    index = np.searchsorted(mass_ratio, z)
    return radii[index], index