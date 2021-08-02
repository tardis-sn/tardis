import copy
import numpy as np
from astropy.coordinates import cartesian_to_spherical

from tardis.energy_input.util import (
    GXPhoton,
    GXPhotonStatus,
    kappa_calculation,
    euler_rodrigues,
    compton_theta_distribution,
    get_random_theta_photon,
    get_random_phi_photon,
    normalize,
    get_perpendicular_vector,
)


def get_compton_angle(energy):
    """
    Computes the compton angle from the Klein-Nishina equation.

    Parameters
    ----------
    energy : float
        Photon energy

    Returns
    -------
    compton_angle : float
        Compton scattering angle
    lost_energy : float
        Energy lost based on angle
    new_energy : float
        Photon energy
    """
    theta_angles, theta_distribution = compton_theta_distribution(energy)

    z = np.random.random()

    # get Compton scattering angle
    compton_angle = np.interp(z, theta_angles, theta_distribution)
    # Energy calculations
    new_energy = energy / (
        1.0 + kappa_calculation(energy) * (1.0 - np.cos(compton_angle))
    )
    lost_energy = energy - new_energy

    return compton_angle, lost_energy, new_energy


def compton_scatter(photon, compton_angle):
    """
    Changes the direction of the gamma-ray by the Compton scattering angle

    Parameters
    ----------
    photon : GXPhoton object
    compton_angle : float

    Returns
    -------
    float64
        Photon theta direction
    float64
        Photon phi direction
    """
    # transform original direction vector to cartesian coordinates
    original_direction = normalize(photon.direction.cartesian_coords)
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

    return theta_final.value + 0.5 * np.pi, phi_final.value


def pair_creation(photon):
    """
    Randomly scatters the input gamma ray
    Sets its energy to 511 KeV
    Creates backwards photon

    Parameters
    ----------
    photon : GXPhoton object

    Returns
    -------
        GXPhoton
            forward photon
        GXPhoton
            backward photon
    """
    direction_theta = get_random_theta_photon()
    direction_phi = get_random_phi_photon()

    photon.energy = 511.0
    photon.direction.theta = direction_theta
    photon.direction.phi = direction_phi

    backward_ray = GXPhoton(
        copy.deepcopy(photon.location),
        copy.deepcopy(photon.direction),
        copy.deepcopy(photon.energy),
        GXPhotonStatus.IN_PROCESS,
        copy.deepcopy(photon.shell),
    )

    backward_ray.direction.phi += np.pi

    if backward_ray.direction.phi > 2 * np.pi:
        backward_ray.direction.phi -= 2 * np.pi

    return photon, backward_ray


def scatter_type(compton_opacity, photoabsorption_opacity, total_opacity):
    """
    Determines the scattering type based on process opacities

    Parameters
    ----------
    compton_opacity : float
    photoabsorption_opacity : float
    total_opacity : float

    Returns
    -------
    status : GXPhotonStatus
        Scattering process the photon encounters

    """
    z = np.random.random()

    if z <= (compton_opacity / total_opacity):
        status = GXPhotonStatus.COMPTON_SCATTER
    elif z <= (compton_opacity + photoabsorption_opacity) / total_opacity:
        status = GXPhotonStatus.PHOTOABSORPTION
    else:
        status = GXPhotonStatus.PAIR_CREATION

    return status
