import astropy.units as u
import tardis.constants as const
import numpy as np

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
    electron_fraction = electron_number_density / number_density

    return electron_fraction * (2 * np.pi * const.e ** 4) / energy * np.log(4 * energy / plasma_frequency)

def cross_section(energy):
    return False

def degradation(energy):
    return False

def heating_fraction(energy, mean_initial_energy, electron_number_density, number_density, lowest_energy, max_energy):
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

    constant = lowest_energy * degradation(lowest_energy) * \
        coulomb_loss_function(lowest_energy, electron_number_density, number_density)

    degradation_energies = np.linspace(lowest_energy, max_energy)
    degradation_spectrum = degradation(degradation_energies) 
    coulomb_loss_spectrum = coulomb_loss_function(degradation_energies, electron_number_density, number_density)

    integral_degradation = np.trapz(degradation_spectrum * coulomb_loss_spectrum)

    number_of_things = 1 * 1

    integral_number = np.trapz(number_of_things)

    return mean_fraction * integral_degradation + mean_fraction * constant + mean_fraction * integral_number

def excitation_fraction(energy, max_energy, number_density, mean_initial_energy):
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
    degradation_spectrum = degradation(degradation_energies)

    cross_sections = cross_section(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * cross_sections)

    return (number_density * energy) / mean_initial_energy * integral_degradation

def ionization_fraction(energy, max_energy, number_density, mean_initial_energy):
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
    degradation_spectrum = degradation(degradation_energies)

    q = q(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * q)

    return (number_density * energy) / mean_initial_energy * integral_degradation
