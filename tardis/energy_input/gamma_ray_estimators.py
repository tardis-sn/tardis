import numpy as np
from numba import njit

from tardis.energy_input.util import (
    ELECTRON_MASS_ENERGY_KEV,
    H_CGS_KEV,
    angle_aberration_gamma,
    doppler_factor_3d,
)
from tardis.opacities.opacities import (
    SIGMA_T,
    compton_opacity_calculation,
    kappa_calculation,
    photoabsorption_opacity_calculation,
)
from tardis.transport.montecarlo import njit_dict_no_parallel


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

    doppler_factor = doppler_factor_3d(
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


@njit(**njit_dict_no_parallel)
def get_average_compton_fraction(energy):
    def f(x, mu):
        return 1.0 / (1.0 + x * (1.0 - mu))

    def cross_section(x, mu):
        return (
            (3.0 * SIGMA_T)
            / (16.0 * np.pi)
            * f(x, mu) ** 2.0
            * (f(x, mu) + 1.0 / f(x, mu) - (1.0 - mu**2))
        )

    x = kappa_calculation(energy)
    mus = np.linspace(-1, 1, 100)
    dmu = mus[1] - mus[0]
    sum = 0
    norm = 0

    for mu in mus:
        sum += cross_section(x, mu) * f(x, mu) * dmu
        norm += cross_section(x, mu) * dmu

    integral = 1.0 - sum / norm

    return 1 - integral


@njit(**njit_dict_no_parallel)
def deposition_estimator_kasen(energy, ejecta_density, iron_group_fraction):
    return get_average_compton_fraction(energy) * compton_opacity_calculation(
        energy, ejecta_density
    ) + photoabsorption_opacity_calculation(
        energy, ejecta_density, iron_group_fraction
    )
