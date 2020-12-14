import astropy.units as u
import tardis.constants as const
import numpy as np

# Eqn 7 set up energy grid, bin-wise integration, multiply S(E) by the energy grid


def coulomb_loss_function(energy, electron_number_density, number_density):
    """Calculates the Coulomb loss function for a given energy,
    electron number density, and atom number density

    Parameters
    ----------
    energy : float
        electron energy
    electron_number_density : float

    number_density : float

    Returns
    -------
    float
        Coulomb loss energy
    """
    plasma_frequency = 56414.6 * np.sqrt(electron_number_density)
    zeta_electron = 2.0 * const.hbar.to_cgs().value * plasma_frequency
    electron_fraction = electron_number_density / number_density

    return (
        electron_fraction
        * (2 * np.pi * const.e ** 4)
        / energy
        * np.log(4 * energy / zeta_electron)
    )


def cross_section(
    energy, initial_electron_energy, oscillator_strength, van_regemorter_fit
):
    # Probably more useful from macro atom
    hydrogen_ionization_potential = 1

    k = initial_electron_energy / 13.60

    return (
        (8 * np.pi)
        / np.sqrt(3)
        * (1 / k ** 2)
        * (hydrogen_ionization_potential / energy)
        * oscillator_strength
        * van_regemorter_fit
        * const.a0 ** 2
    )


def collisional_cross_section(energy):
    # Younger 1981 apparently
    # But also used in macro atom
    return False


def electron_spectrum(energy):
    """number of electrons at energy E

    Parameters
    ----------
    energy : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return False


def heating_fraction(
    energy,
    mean_initial_energy,
    electron_number_density,
    number_density,
    lowest_energy,
    max_energy,
):
    """Calculates the heating fraction of electrons

    Parameters
    ----------
    energy : float

    mean_initial_energy : float

    electron_number_density : float

    number_density : float

    lowest_energy : float

    max_energy : float

    Returns
    -------
    float
        fraction of energy going into heating
    """
    mean_fraction = 1 / mean_initial_energy

    constant = (
        lowest_energy
        * electron_spectrum(lowest_energy)
        * coulomb_loss_function(
            lowest_energy, electron_number_density, number_density
        )
    )

    degradation_energies = np.linspace(lowest_energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)
    coulomb_loss_spectrum = coulomb_loss_function(
        degradation_energies, electron_number_density, number_density
    )

    integral_degradation = np.trapz(
        degradation_spectrum * coulomb_loss_spectrum
    )

    number_of_things = 1 * 1

    integral_number = np.trapz(number_of_things)

    return (
        mean_fraction * integral_degradation
        + mean_fraction * constant
        + mean_fraction * integral_number
    )


def excitation_fraction(
    energy, max_energy, number_density, mean_initial_energy
):
    """Excitation fraction of electron energy

    Parameters
    ----------
    energy : float
        [description]
    max_energy : float
        [description]
    number_density : float
        [description]
    mean_initial_energy : float
        [description]

    Returns
    -------
    float
        fraction of energy going into excitation
    """
    degradation_energies = np.linspace(energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)

    cross_sections = cross_section(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * cross_sections)

    return (
        (number_density * energy) / mean_initial_energy * integral_degradation
    )


def ionization_fraction(
    energy, max_energy, number_density, mean_initial_energy
):
    """Ionization fraction of electron energy

    Parameters
    ----------
    energy : float
        [description]
    max_energy : float
        [description]
    number_density : float
        [description]
    mean_initial_energy : float
        [description]

    Returns
    -------
    float
        Ionization fraction of electron energy
    """
    degradation_energies = np.linspace(energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)

    q = collisional_cross_section(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * q)

    return (
        (number_density * energy) / mean_initial_energy * integral_degradation
    )
