import numpy as np
import pandas as pd

# Need to scale deposited energy to rate of decay reaction to get energy per second in steady state
# Assume photon path is small compared to dynamical effects


def read_nuclear_dataframe(path):
    """Reads HDF5 output from tardis nuclear

    Parameters
    ----------
    path : str
        Path to HDF5 file

    Returns
    -------
    Pandas dataframe
        Dataframe of decay radiation properties
    """
    return pd.read_hdf(path, key="decay_radiation")


def get_type_property(nuclear_df, type_of_radiation, property):
    """Queries a dataframe of decay radiation properties

    Parameters
    ----------
    nuclear_df : Pandas dataframe
        Dataframe of nuclear decay properties
    type_of_radiation : str
        The type of radiation to get properties from
    property : str
        Property to return

    Returns
    -------
    One-dimensional Numpy Array, dtype float
        Array of values matching the queried type and property
    """
    return nuclear_df.query("type==" + type_of_radiation)[property].values


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
    norm_intensity = intensity / np.sum(intensity)
    cdf = np.zeros_like(norm_intensity)

    for index, i in enumerate(norm_intensity):
        cdf[index] = cdf[index - 1] + i

    energy.sort()

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
    z = np.random.random()

    index = np.searchsorted(cdf, z)

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
    intensity = get_type_property(nuclear_data, source, "intensity")
    energy = get_type_property(nuclear_data, source, "energy")
    energy_sorted, cdf = create_energy_cdf(energy, intensity)

    return energy_sorted, cdf


def intensity_ratio(nuclear_data, source_1, source_2):
    """Determined the ratio of intensities between two
    sources of decay radiation

    Parameters
    ----------
    nuclear_data : Pandas dataframe
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
    intensity_1 = get_type_property(nuclear_data, source_1, "intensity")
    intensity_2 = get_type_property(nuclear_data, source_2, "intensity")
    total_intensity = np.sum(intensity_1) + np.sum(intensity_2)
    return (
        np.sum(intensity_1) / total_intensity,
        np.sum(intensity_2) / total_intensity,
    )
