import copy
from typing_extensions import final
import numpy as np
from numba import njit

from tardis.energy_input.util import (
    doppler_gamma,
    kappa_calculation,
    euler_rodrigues,
    compton_theta_distribution,
    get_random_theta_photon,
    get_random_phi_photon,
    normalize_vector,
    get_perpendicular_vector,
    cartesian_to_spherical,
    angle_aberration_gamma,
    ELECTRON_MASS_ENERGY_KEV,
)
from tardis.energy_input.GXPhoton import GXPhotonStatus


@njit
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
    compton_angle = np.interp(z, theta_distribution, theta_angles)
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
        Photon energy
    float64
        Photon theta direction
    float64
        Photon phi direction
    """
    # transform original direction vector to cartesian coordinates
    original_direction = normalize_vector(
        np.array(photon.direction.cartesian_coords)
    )

    # get comoving frame direction
    comov_direction = angle_aberration_gamma(
        photon.direction.vector, photon.location.r
    )

    # compute an arbitrary perpendicular vector to the original direction
    orthogonal_vector = get_perpendicular_vector(original_direction)
    # determine a random vector with compton_angle to the original direction
    # TODO is this correctly the comoving frame direction? Or should this be rest frame?
    new_vector = normalize_vector(
        np.dot(
            euler_rodrigues(compton_angle, orthogonal_vector),
            comov_direction,
        )
    )

    comoving_energy = photon.energy * doppler_gamma(
        photon.direction.vector, photon.location.r
    )

    # draw a random angle from [0,2pi]
    phi = 2.0 * np.pi * np.random.random()
    # rotate the vector with compton_angle around the original direction
    # TODO is this correctly the comoving frame direction? Or should this be rest frame?
    final_compton_scattered_vector = normalize_vector(
        np.dot(euler_rodrigues(phi, comov_direction), new_vector)
    )
    # transform the cartesian coordinates to spherical coordinates
    final_comov_direction = np.array(
        cartesian_to_spherical(
            final_compton_scattered_vector[0],
            final_compton_scattered_vector[1],
            final_compton_scattered_vector[2],
        )
    )

    # Calculate the angle aberration of the final direction
    final_direction = angle_aberration_gamma(
        final_comov_direction, photon.location.r
    )

    # Transform the energy back to the lab frame
    energy = comoving_energy / doppler_gamma(final_direction, photon.location.r)
    return energy, final_direction[1], final_direction[2]


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

    # Calculate aberration of the random angle for the rest frame
    final_direction = angle_aberration_gamma(
        np.array([photon.direction.r, direction_theta, direction_phi]),
        photon.location.r,
    )

    # TODO Is this correct? Scaling photon energy by the inverse doppler factor
    photon.energy = ELECTRON_MASS_ENERGY_KEV / doppler_gamma(
        final_direction, photon.location.r
    )
    photon.direction.theta = final_direction[1]
    photon.direction.phi = final_direction[2]

    backward_ray = copy.deepcopy(photon)
    backward_ray.status = GXPhotonStatus.IN_PROCESS

    # TODO Is this correct? Should this have aberration
    backward_ray.direction.phi += np.pi

    if backward_ray.direction.phi > 2 * np.pi:
        backward_ray.direction.phi -= 2 * np.pi

    return photon, backward_ray


@njit
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
