import numpy as np
import pandas as pd
from nuclear.ejecta import Ejecta


def decay_nuclides(shell_mass, initial_composition, epoch):
    decay_model = Ejecta(shell_mass, initial_composition)

    new_fractions = decay_model.decay(epoch)
    return new_fractions


def create_energy_cdf(energy, intensity):
    """Creates a CDF of given intensities

    Parameters
    ----------
    energy :  One-dimensional Numpy Array, dtype float
        Array of energies
    intensity :  One-dimensional Numpy Array, dtype float
        Array of intensities

    Returns
    -------
    One-dimensional Numpy Array, dtype float
        Sorted energy array
    One-dimensional Numpy Array, dtype float
        CDF where each index corresponds to the energy in
        the sorted array
    """
    energy.sort()
    sorted_indices = np.argsort(energy)
    sorted_intensity = intensity[sorted_indices]
    norm_intensity = sorted_intensity / np.sum(sorted_intensity)
    cdf = np.cumsum(norm_intensity)

    return energy, cdf


def sample_energy_distribution(energy_sorted, cdf):
    """Randomly samples a CDF of energies

    Parameters
    ----------
    energy_sorted : One-dimensional Numpy Array, dtype float
        Sorted energy array
    cdf : One-dimensional Numpy Array, dtype float
        CDF

    Returns
    -------
    float
        Sampled energy
    """
    index = np.searchsorted(cdf, np.random.random())

    return energy_sorted[index]


def setup_input_energy(nuclear_data, source):
    """Sets up energy distribution and CDF for a
    source of decay radiation.

    Parameters
    ----------
    nuclear_data : Pandas dataframe
        Dataframe of nuclear decay properties
    source : str
        Type of decay radiation

    Returns
    -------
    One-dimensional Numpy Array, dtype float
        Sorted energy array
    One-dimensional Numpy Array, dtype float
        CDF where each index corresponds to the energy in
        the sorted array
    """
    intensity = nuclear_data.query("type==" + source)["intensity"].values
    energy = nuclear_data.query("type==" + source)["energy"].values
    energy_sorted, cdf = create_energy_cdf(energy, intensity)

    return energy_sorted, cdf


def intensity_ratio(nuclear_data, source_1, source_2):
    """Determined the ratio of intensities between two
    sources of decay radiation

    Parameters
    ----------
    nuclear_data : pandas.Dataframe
        Dataframe of nuclear decay properties
    source_1 : str
        Type of decay radiation to compare
    source_2 : str
        Type of decay radiation to compare

    Returns
    -------
    float
        Fractional intensity of source_1
    float
        Fractional intensity of source_2
    """
    intensity_1 = nuclear_data.query("type==" + source_1)["intensity"].values
    intensity_2 = nuclear_data.query("type==" + source_2)["intensity"].values
    total_intensity = np.sum(intensity_1) + np.sum(intensity_2)
    scale_factor = total_intensity / 100
    print("Scale factor: ", scale_factor)
    return (
        np.sum(intensity_1) / total_intensity,
        np.sum(intensity_2) / total_intensity,
        scale_factor,
    )
