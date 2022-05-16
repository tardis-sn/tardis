import numpy as np

from tardis.energy_input.util import (
    angle_aberration_gamma,
    doppler_gamma,
    H_CGS_KEV,
    ELECTRON_MASS_ENERGY_KEV,
)
from tardis.energy_input.calculate_opacity import (
    pair_creation_opacity_calculation,
    SIGMA_T,
)


def compton_emissivity_estimator(packet, distance):
    """Compton scattering emissivity estimator for integral
    calculations

    Parameters
    ----------
    packet : GXPacket
        Packet that needs its emissivity calculated
    distance : float64
        Distance packet has travelled

    Returns
    -------
    float64, int
        Unnormalized emissivity estimator, line index
    """

    cmf_direction = angle_aberration_gamma(
        packet.get_direction_vector(), packet.location_r
    )

    cmf_angle = np.dot(cmf_direction, [1, 0, 0])

    frequency_factor = (
        1
        + H_CGS_KEV * packet.nu_cmf / ELECTRON_MASS_ENERGY_KEV
        + (1.0 - cmf_angle)
    )

    line_index = GET_NEAREST_LINE_REDWARD_FUNCTION(
        packet.nu_cmf / frequency_factor
    )

    partial_cross_section = (
        3.0
        / 16.0
        / np.pi
        * SIGMA_T
        / frequency_factor
        / frequency_factor
        * (frequency_factor + (1.0 / frequency_factor) + cmf_angle**2.0 - 1.0)
    )

    doppler_factor = doppler_gamma(
        packet.get_direction_vector(), packet.location_r
    )

    emissivity = (
        packet.energy_rf
        * partial_cross_section
        * distance
        * doppler_factor**2.0
        / frequency_factor
    )

    return emissivity, line_index


def pair_creation_estimator(packet, pair_creation_opacity, distance):
    """Calculates the emissivity for pair creation gamma-rays

    Parameters
    ----------
    packet : GXPacket
        Packet that needs its emissivity calculated
    pair_creation_opacity : float64
        Opacity of the pair creation process
    distance : float64
        Distance packet has travelled

    Returns
    -------
    float64
        Emissivity estimator
    """
    normalized_energy = (
        2 * ELECTRON_MASS_ENERGY_KEV / (H_CGS_KEV * packet.nu_cmf)
    )

    emissivity = (
        pair_creation_opacity * normalized_energy * packet.energy_rf * distance
    )

    return emissivity
