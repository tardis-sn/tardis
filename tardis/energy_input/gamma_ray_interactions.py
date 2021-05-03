import numpy as np
from astropy.coordinates import cartesian_to_spherical
from astropy.coordinates import spherical_to_cartesian
from tardis.energy_input.util import (
    kappa_calculation,
    klein_nishina,
    euler_rodrigues,
    compton_theta_distribution,
    SphericalVector,
    get_random_theta_gamma_ray,
    get_random_phi_gamma_ray,
)

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu


def get_compton_angle(gamma_ray):
    """
    Computes the compton angle from the Klein-Nishina equation.

    Parameters
    ----------
    gamma_ray : GammaRay object

    Returns
    -------
    compton_angle : dtype float
    lost_energy : dtype float
    """

    theta_angles, theta_distribution = compton_theta_distribution(
        gamma_ray.energy
    )

    z = np.random.random()

    # get Compton scattering angle
    compton_angle = theta_angles[np.searchsorted(theta_distribution, z)]
    # Energy calculations
    new_energy = gamma_ray.energy / (
        1.0
        + kappa_calculation(gamma_ray.energy) * (1.0 - np.cos(compton_angle))
    )
    lost_energy = gamma_ray.energy - new_energy
    gamma_ray.energy = new_energy

    return compton_angle, lost_energy


def compton_scatter(gamma_ray, compton_angle):
    """
    Changes the direction of the gamma-ray by the Compton scattering angle

    Parameters
    ----------
    gamma_ray : GammaRay object
    compton_angle : dtype float

    Returns
    -------

    """
    # transform original direction vector to cartesian coordinates
    original_direction = normalize(gamma_ray.direction.get_cartesian_coords)
    # compute an arbitrary perpendicular vector to the original direction
    orthogonal_vector = get_perpendicular_vector(original_direction)
    # determine a random vector with compton_angle to the original direction
    new_vector = normalize(
        np.dot(
            euler_rodrigues(compton_angle, orthogonal_vector),
            original_direction,
        )
    )

    # draw a random angle from [0,2pi]
    phi = 2.0 * np.pi * np.random.random()
    # rotate the vector with compton_angle around the original direction
    final_compton_scattered_vector = normalize(
        np.dot(euler_rodrigues(phi, original_direction), new_vector)
    )
    # transform the cartesian coordinates to spherical coordinates
    r_final, theta_final, phi_final = cartesian_to_spherical(
        final_compton_scattered_vector[0],
        final_compton_scattered_vector[1],
        final_compton_scattered_vector[2],
    )
    # 0.5*np.pi added because of the definition of theta
    # in astropy.coordinates.cartesian_to_spherical
    gamma_ray.direction.theta = theta_final.value + 0.5 * np.pi
    gamma_ray.direction.phi = phi_final.value


def normalize(vector):
    """
    Normalizes a vector in cartesian coordinates

    Parameters
    ----------
    vector : One-dimensional Numpy Array, dtype float

    Returns
    -------
    normalized_vector : One-dimensional Numpy Array, dtype float
    """
    normalized_vector = vector / np.sqrt(
        vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
    )
    return normalized_vector


def get_perpendicular_vector(original_direction):
    """
    Computes a vector which is perpendicular to the input vector

    Parameters
    ----------
    original_direction : SphericalVector object

    Returns
    -------

    """
    # draw random angles
    theta = get_random_theta_gamma_ray()
    phi = get_random_phi_gamma_ray()
    # transform random angles to cartesian coordinates
    # 0.5*np.pi subtracted because of the definition of theta
    # in astropy.coordinates.cartesian_to_spherical
    random_vector = spherical_to_cartesian(1, theta - 0.5 * np.pi, phi)
    perpendicular_vector = np.cross(original_direction, random_vector)
    perpendicular_vector = normalize(perpendicular_vector)
    return perpendicular_vector


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
    direction_theta = get_random_theta_gamma_ray()
    direction_phi = get_random_phi_gamma_ray()

    gamma_ray.energy = 511.0
    gamma_ray.direction.theta = direction_theta
    gamma_ray.direction.phi = direction_phi


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
    compton_angle : dtype float

    """
    z = np.random.random()

    ejecta_energy_gain = 0.0
    compton_angle = 0.0

    if z <= (compton_opacity / total_opacity):
        gamma_ray.status = "ComptonScatter"
        # not happy about computing the compton angle twice.
        # should be restructured
        compton_angle, ejecta_energy_gain = get_compton_angle(gamma_ray)
    elif z <= (compton_opacity + photoabsorption_opacity) / total_opacity:
        gamma_ray.status = "PhotoAbsorbed"
        ejecta_energy_gain = gamma_ray.energy
    else:
        gamma_ray.status = "PairCreated"
        ejecta_energy_gain = gamma_ray.energy - (2.0 * 511.0)

    return ejecta_energy_gain, compton_angle
