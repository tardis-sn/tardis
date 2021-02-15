import numpy as np
import pandas as pd

# Need to scale deposited energy to rate of decay reaction to get energy per second in steady state
# Assume photon path is small compared to dynamical effects


def read_nuclear_dataframe(path):
    return pd.read_hdf(path, key="decay_radiation")


def get_type_property(nuclear_df, type_of_radiation, property):
    return nuclear_df.query("type==" + type_of_radiation)[property].values


def create_energy_cdf(energy, intensity):
    norm_intensity = intensity / np.sum(intensity)
    cdf = np.zeros_like(norm_intensity)

    for index, i in enumerate(norm_intensity):
        cdf[index] = cdf[index - 1] + i

    energy.sort()

    return energy, cdf


def sample_energy_distribution(energy_sorted, cdf):

    z = np.random.random()

    index = np.searchsorted(cdf, z)

    return energy_sorted[index]


def setup_input_energy(nuclear_data, source):
    intensity = get_type_property(nuclear_data, source, "intensity")
    energy = get_type_property(nuclear_data, source, "energy")
    energy_sorted, cdf = create_energy_cdf(energy, intensity)

    return energy_sorted, cdf


def intensity_ratio(nuclear_data, source_1, source_2):
    intensity_1 = get_type_property(nuclear_data, source_1, "intensity")
    intensity_2 = get_type_property(nuclear_data, source_2, "intensity")
    total_intensity = np.sum(intensity_1) + np.sum(intensity_2)
    return (
        np.sum(intensity_1) / total_intensity,
        np.sum(intensity_2) / total_intensity,
    )
