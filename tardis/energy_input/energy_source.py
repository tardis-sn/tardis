import numpy as np
import pandas as pd

def read_nuclear_dataframe(path):
    return pd.read_hdf(path, key="decay_radiation")

def get_gamma_rays_property(nuclear_df, property):
    return nuclear_df.query("type=='gamma_rays'")[property].values

def create_energy_cdf(energy, intensity):
    norm_intensity = intensity / np.sum(intensity)
    cdf = np.zeros_like(norm_intensity)

    for index, i in enumerate(norm_intensity):
        cdf[index] = cdf[index-1] + i

    energy.sort()

    return energy, cdf

def sample_energy_distribution(energy_sorted, cdf):

    z = np.random.random()

    index = np.searchsorted(cdf, z)

    return energy_sorted[index]

def setup_gamma_ray_energy(nuclear_data):
    intensity = get_gamma_rays_property(nuclear_data, "intensity")
    energy = get_gamma_rays_property(nuclear_data, "energy")
    energy_sorted, cdf = create_energy_cdf(energy, intensity)

    return energy_sorted, cdf