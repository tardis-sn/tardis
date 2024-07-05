import numpy as np
from numba import njit

from tardis.energy_input.GXPacket import GXPacketStatus
from tardis.energy_input.util import (
    ELECTRON_MASS_ENERGY_KEV,
    H_CGS_KEV,
    angle_aberration_gamma,
    compton_theta_distribution,
    doppler_factor_3d,
    euler_rodrigues,
    get_perpendicular_vector,
    get_random_unit_vector,
)
from tardis.opacities.opacities import (
    compton_opacity_partial,
    kappa_calculation,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


@njit(**njit_dict_no_parallel)
def get_compton_angle(energy):
    """
    Computes the compton angle from the Klein-Nishina equation.
    Computes the lost energy due to this angle

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


@njit(**njit_dict_no_parallel)
def get_compton_fraction(energy):
    """
    Computes the compton angle from the Klein-Nishina equation.
    Determines the probability of absorption from this angle.

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


@njit(**njit_dict_no_parallel)
def get_compton_fraction_artis(energy):
    """Gets the Compton scattering/absorption fraction
    and angle following the scheme in ARTIS

    Parameters
    ----------
    energy : float
        Energy of the gamma-ray

    Returns
    -------
    float
        Scattering angle
    float
        Compton scattering fraction
    """
    energy_norm = kappa_calculation(energy)

    fraction_max = 1.0 + 2.0 * energy_norm
    fraction_min = 1.0

    normalization = np.random.random() * compton_opacity_partial(
        energy_norm, fraction_max
    )

    epsilon = 1.0e20
    count = 0

    while epsilon > 1.0e-4:
        fraction_try = (fraction_max + fraction_min) / 2.0
        sigma_try = compton_opacity_partial(energy_norm, fraction_try)

        if sigma_try > normalization:
            fraction_max = fraction_try
            epsilon = (sigma_try - normalization) / normalization
        else:
            fraction_min = fraction_try
            epsilon = (normalization - sigma_try) / normalization

        count += 1
        if count > 1000:
            print("Error, failure to get a Compton fraction")
            break

    angle = np.arccos(1.0 - ((fraction_try - 1) / energy_norm))

    return angle, fraction_try


@njit(**njit_dict_no_parallel)
def get_compton_fraction_urilight(energy):
    """Gets the Compton scattering/absorption fraction
    and angle following the scheme in Urilight

    Parameters
    ----------
    energy : float
        Energy of the gamma-ray

    Returns
    -------
    float
        Scattering angle
    float
        Compton scattering fraction
    """
    E0 = kappa_calculation(energy)

    x0 = 1.0 / (1.0 + 2.0 * E0)

    accept = False
    while not accept:
        z = np.random.random(3)
        alpha1 = np.log(1.0 / x0)
        alpha2 = (1.0 - x0**2.0) / 2.0
        if z[1] < alpha1 / (alpha1 + alpha2):
            x = x0 ** z[2]
        else:
            x = np.sqrt(x0**2.0 + (1.0 - x0**2.0) * z[2])

        f = (1.0 - x) / x / E0
        sin2t = f * (2.0 - f)
        ge = 1.0 - x / (1 + x**2.0) * sin2t
        if ge > z[3]:
            accept = True

    cost = 1.0 - f

    return np.arccos(cost), x


@njit(**njit_dict_no_parallel)
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
    comov_direction = angle_aberration_gamma(
        photon.direction, photon.location, photon.time_current
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

    # Calculate the angle aberration of the final direction
    final_direction = angle_aberration_gamma(
        final_compton_scattered_vector, photon.location, photon.time_current
    )

    return final_direction


@njit(**njit_dict_no_parallel)
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

    new_direction = get_random_unit_vector()

    # Calculate aberration of the random angle for the rest frame
    final_direction = angle_aberration_gamma(
        new_direction, packet.location, -1 * packet.time_current
    )

    packet.direction = final_direction

    doppler_factor = doppler_factor_3d(
        packet.direction, packet.location, packet.time_current
    )

    packet.nu_cmf = ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    packet.nu_rf = packet.nu_cmf / doppler_factor
    packet.energy_rf = packet.energy_cmf / doppler_factor

    return packet


@njit(**njit_dict_no_parallel)
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
    status : GXPacketStatus
        Scattering process the photon encounters

    """
    z = np.random.random()

    if z <= (compton_opacity / total_opacity):
        status = GXPacketStatus.COMPTON_SCATTER
    elif z <= (compton_opacity + photoabsorption_opacity) / total_opacity:
        status = GXPacketStatus.PHOTOABSORPTION
    else:
        status = GXPacketStatus.PAIR_CREATION

    return status
