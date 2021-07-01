import copy
import numpy as np
from astropy.coordinates import cartesian_to_spherical
from astropy.coordinates import spherical_to_cartesian

from tardis.energy_input.base import GXPacketStatus, GXPacket
from tardis.energy_input.util import (
    kappa_calculation,
    euler_rodrigues,
    compton_theta_distribution,
    get_random_theta_gamma_ray,
    get_random_phi_gamma_ray,
    normalize,
    get_perpendicular_vector,
)

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu


def get_compton_angle(gxpacket):
    """
    Computes the compton angle from the Klein-Nishina equation.

    Parameters
    ----------
    gxpacket : GXPacket object

    Returns
    -------
    compton_angle : float
    lost_energy : float
    """

    theta_angles, theta_distribution = compton_theta_distribution(
        gxpacket.energy
    )

    z = np.random.random()

    # get Compton scattering angle
    compton_angle = theta_angles[np.searchsorted(theta_distribution, z)]
    # Energy calculations
    new_energy = gxpacket.energy / (
        1.0 + kappa_calculation(gxpacket.energy) * (1.0 - np.cos(compton_angle))
    )
    lost_energy = gxpacket.energy - new_energy
    gxpacket.energy = new_energy

    return compton_angle, lost_energy


def compton_scatter(gxpacket, compton_angle):
    """
    Changes the direction of the gamma-ray by the Compton scattering angle

    Parameters
    ----------
    gxpacket : GXPacket object
    compton_angle : float

    Returns
    -------
    float64
        Packet theta direction
    float64
        Packet phi direction
    """
    # transform original direction vector to cartesian coordinates
    original_direction = normalize(gxpacket.direction.cartesian_coords)
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


def pair_creation(gxpacket):
    """
    Randomly scatters the input gamma ray
    Sets its energy to 511 KeV
    Creates backwards packet

    Parameters
    ----------
    gxpacket : GXPacket object

    Returns
    -------
        GXPacket
            forward packet
        GXPacket
            backward packet
    """
    direction_theta = get_random_theta_gamma_ray()
    direction_phi = get_random_phi_gamma_ray()

    gxpacket.energy = 511.0
    gxpacket.direction.theta = direction_theta
    gxpacket.direction.phi = direction_phi

    backward_ray = GXPacket(
        copy.deepcopy(gxpacket.location),
        copy.deepcopy(gxpacket.direction),
        copy.deepcopy(gxpacket.energy),
        GXPacketStatus.IN_PROCESS,
        copy.deepcopy(gxpacket.shell),
    )

    backward_ray.direction.phi += np.pi

    if backward_ray.direction.phi > 2 * np.pi:
        backward_ray.direction.phi -= 2 * np.pi

    return gxpacket, backward_ray


def scatter_type(
    gxpacket, compton_opacity, photoabsorption_opacity, total_opacity
):
    """
    Determines the scattering type based on process opacities

    Parameters
    ----------
    gxpacket : GXPacket object
    compton_opacity : float
    photoabsorption_opacity : float
    total_opacity : float

    Returns
    -------
    ejecta_energy_gain : float
    compton_angle : float

    """
    z = np.random.random()

    ejecta_energy_gain = 0.0
    compton_angle = 0.0

    if z <= (compton_opacity / total_opacity):
        gxpacket.status = GXPacketStatus.COMPTON_SCATTER
        compton_angle, ejecta_energy_gain = get_compton_angle(gxpacket)
    elif z <= (compton_opacity + photoabsorption_opacity) / total_opacity:
        gxpacket.status = GXPacketStatus.PHOTOABSORPTION
        ejecta_energy_gain = gxpacket.energy
    else:
        gxpacket.status = GXPacketStatus.PAIR_CREATION
        ejecta_energy_gain = gxpacket.energy - (2.0 * 511.0)

    return ejecta_energy_gain, compton_angle
