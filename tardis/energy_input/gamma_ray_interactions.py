import numpy as np
from numba import njit

from tardis.energy_input.util import (
    kappa_calculation,
    euler_rodrigues,
    compton_theta_distribution,
    get_random_theta_photon,
    get_random_phi_photon,
    get_perpendicular_vector,
    cartesian_to_spherical,
    angle_aberration_gamma,
    spherical_to_cartesian,
    doppler_gamma,
    ELECTRON_MASS_ENERGY_KEV,
    H_CGS_KEV,
)
from tardis.energy_input.GXPhoton import GXPhotonStatus, GXPhoton
from tardis.energy_input.GXPacket import GXPacketStatus


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


@njit
def get_compton_fraction(energy):
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
    compton_fraction : float
        Fraction of energy lost
    """
    theta_angles, theta_distribution = compton_theta_distribution(energy)

    z = np.random.random()

    # get Compton scattering angle
    compton_angle = np.interp(z, theta_distribution, theta_angles)
    # Energy calculations
    fraction = 1 / (
        1.0 + kappa_calculation(energy) * (1.0 - np.cos(compton_angle))
    )

    return compton_angle, fraction


@njit
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

    # get comoving frame direction
    abb_array = angle_aberration_gamma(
        photon.get_direction_vector(), photon.location_r
    )
    comov_direction = np.array(
        spherical_to_cartesian(1, abb_array[1], abb_array[2])
    )

    # compute an arbitrary perpendicular vector to the comoving direction
    orthogonal_vector = get_perpendicular_vector(comov_direction)
    # determine a random vector with compton_angle to the comoving direction
    new_vector = np.dot(
        euler_rodrigues(compton_angle, orthogonal_vector),
        comov_direction,
    )

    # draw a random angle from [0,2pi]
    phi = 2.0 * np.pi * np.random.random()
    # rotate the vector with compton_angle around the comoving direction
    final_compton_scattered_vector = np.dot(
        euler_rodrigues(phi, comov_direction), new_vector
    )

    norm_phi = np.dot(
        final_compton_scattered_vector, final_compton_scattered_vector
    )

    norm_theta = np.dot(final_compton_scattered_vector, comov_direction)

    assert (
        np.abs(norm_phi - 1) < 1e-8
    ), "Error, norm of Compton scatter vector is not 1!"

    assert (
        np.abs(norm_theta - np.cos(compton_angle)) < 1e-8
    ), "Error, difference between new vector angle and Compton angle is more than 0!"

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
        final_comov_direction, photon.location_r
    )

    return final_direction[1], final_direction[2]


@njit
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
        np.array([1.0, direction_theta, direction_phi]),
        photon.location_r,
    )

    photon.energy = ELECTRON_MASS_ENERGY_KEV
    photon.direction_theta = final_direction[1]
    photon.direction_phi = final_direction[2]

    backward_ray = GXPhoton(
        photon.location_r,
        photon.location_theta,
        photon.location_phi,
        final_direction[1],
        final_direction[2],
        photon.energy,
        GXPhotonStatus.IN_PROCESS,
        photon.shell,
        photon.activity,
    )

    # TODO Is this correct? Should this have aberration
    backward_ray.direction_phi += np.pi
    backward_ray.tau = photon.tau
    backward_ray.time_created = photon.time_created
    backward_ray.time_current = photon.time_current

    if backward_ray.direction_phi > 2 * np.pi:
        backward_ray.direction_phi -= 2 * np.pi

    return photon, backward_ray


@njit
def pair_creation_packet(packet):
    """
    Pair creation randomly scatters the packet
    or destroys it, based on the frequency

    Parameters
    ----------
    packet : GXPacket
        incoming packet

    Returns
    -------
    GXPacket
        outgoing packet
    """

    probability_gamma = (
        2 * ELECTRON_MASS_ENERGY_KEV / (H_CGS_KEV * packet.nu_cmf)
    )

    if np.random.random() > probability_gamma:
        packet.status = GXPacketStatus.PHOTOABSORPTION
        return packet

    direction_theta = get_random_theta_photon()
    direction_phi = get_random_phi_photon()

    # Calculate aberration of the random angle for the rest frame
    final_direction = angle_aberration_gamma(
        np.array([1.0, direction_theta, direction_phi]),
        packet.location_r,
    )

    packet.direction_theta = final_direction[1]
    packet.direction_phi = final_direction[2]

    doppler_factor = doppler_gamma(
        packet.get_direction_vector(),
        packet.location_r,
    )

    packet.nu_cmf = ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    packet.nu_rf = packet.nu_cmf / doppler_factor
    packet.energy_rf = packet.energy_cmf / doppler_factor

    return packet


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
