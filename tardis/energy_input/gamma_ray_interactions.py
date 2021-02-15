import numpy as np

from tardis.energy_input.util import (
    kappa_calculation,
    klein_nishina,
    euler_rodrigues,
    compton_theta_distribution,
    SphericalVector,
    get_random_mu_gamma_ray,
    get_random_phi_gamma_ray,
)

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu


def compton_scatter(gamma_ray):
    """
    Gets the Compton scattering angle for the input gamma ray
    Calculates the lost energy due to the scattering

    Parameters
    ----------
    gamma_ray : GammaRay object

    Returns
    -------
    lost_energy : dtype float

    """
    theta_angles, theta_distribution = compton_theta_distribution(
        gamma_ray.energy
    )

    z = np.random.random()

    # sample new random theta direction
    random_vector = SphericalVector(
        1.0, get_random_mu_gamma_ray(), get_random_phi_gamma_ray()
    )

    perpendicular_vector_x = (
        gamma_ray.direction.y * random_vector.z
        - gamma_ray.direction.z * random_vector.y
    )
    perpendicular_vector_y = (
        gamma_ray.direction.z * random_vector.x
        - gamma_ray.direction.x * random_vector.z
    )
    perpendicular_vector_z = (
        gamma_ray.direction.x * random_vector.y
        - gamma_ray.direction.y * random_vector.x
    )

    # perpendicular_vector = np.array(
    #    [perpendicular_vector_x, perpendicular_vector_y, perpendicular_vector_z]
    # )

    perpendicular_vector = SphericalVector(
        1.0, perpendicular_vector_z, np.arccos(perpendicular_vector_x)
    )

    # get Compton scattering angle
    compton_angle = theta_angles[np.searchsorted(theta_distribution, z)]

    # rotate to match
    rotation_matrix = euler_rodrigues(compton_angle, perpendicular_vector)
    resulting_direction = np.dot(
        rotation_matrix,
        np.array(
            (
                gamma_ray.direction.x,
                gamma_ray.direction.y,
                gamma_ray.direction.z,
            )
        ),
    )

    phi = get_random_phi_gamma_ray()
    rotation_matrix_phi = euler_rodrigues(phi, gamma_ray.direction)
    final_direction = np.dot(rotation_matrix_phi, resulting_direction)

    gamma_ray.direction.r = 1.0
    gamma_ray.direction.mu = final_direction[2]
    gamma_ray.direction.phi = np.arccos(final_direction[1])

    if gamma_ray.direction.mu < -1.0:
        gamma_ray.direction.mu += 2

    if gamma_ray.direction.phi < 0.0:
        gamma_ray.direction.phi += 2 * np.pi

    # Energy calculations
    new_energy = gamma_ray.energy / (
        1.0
        + kappa_calculation(gamma_ray.energy) * (1.0 - np.cos(compton_angle))
    )
    lost_energy = gamma_ray.energy - new_energy
    gamma_ray.energy = new_energy

    return lost_energy


def pair_creation(gamma_ray):
    """
    Randomly scatters the input gamma ray
    Sets its energy to 511 KeV

    Parameters
    ----------
    gamma_ray : GammaRay object

    Returns
    -------

    """
    direction_mu = get_random_mu_gamma_ray()

    gamma_ray.energy = 511.0
    gamma_ray.direction.mu = direction_mu


def scatter_type(
    gamma_ray, compton_opacity, photoabsorption_opacity, total_opacity
):
    """
    Determines the scattering type based on process opacities

    Parameters
    ----------
    gamma_ray : GammaRay object
    compton_opacity : dtype float
    photoabsorption_opacity : dtype float
    total_opacity : dtype float

    Returns
    -------
    ejecta_energy_gain : dtype float
    pair_created : dtype bool

    """
    z = np.random.random()

    ejecta_energy_gain = 0.0

    if z <= (compton_opacity / total_opacity):
        gamma_ray.status = "ComptonScatter"
        ejecta_energy_gain = compton_scatter(gamma_ray)
    elif z <= (compton_opacity + photoabsorption_opacity) / total_opacity:
        gamma_ray.status = "PhotoAbsorbed"
        ejecta_energy_gain = gamma_ray.energy
    else:
        gamma_ray.status = "PairCreated"
        pair_creation(gamma_ray)

    return ejecta_energy_gain