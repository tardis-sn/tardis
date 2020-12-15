import numpy as np
import pandas as pd


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


def setup_gamma_ray_energy(nuclear_data):
    intensity = get_type_property(nuclear_data, "'gamma_rays'", "intensity")
    intensity = np.append(
        intensity,
        [np.sum(get_type_property(nuclear_data, "'e+'", "intensity"))],
    )
    energy = get_type_property(nuclear_data, "'gamma_rays'", "energy")
    energy = np.append(energy, [511.0])
    energy_sorted, cdf = create_energy_cdf(energy, intensity)

    return energy_sorted, cdf