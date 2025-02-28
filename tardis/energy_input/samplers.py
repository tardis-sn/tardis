import astropy.constants as const
from scipy.interpolate import interp1d
import astropy.units as u
import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel


@njit(**njit_dict_no_parallel)
def sample_mass(masses, inner_radius, outer_radius):
    """Samples location weighted by mass

    Parameters
    ----------
    masses : array
        Shell masses
    inner_radius : array
        Inner radii
    outer_radius : array
        Outer radii

    Returns
    -------
    float
        Sampled radius
    int
        Sampled shell index
    """
    norm_mass = masses / np.sum(masses)
    cdf = np.cumsum(norm_mass)
    shell = np.searchsorted(cdf, np.random.random())

    z = np.random.random()
    radius = (
        z * inner_radius[shell] ** 3.0 + (1.0 - z) * outer_radius[shell] ** 3.0
    ) ** (1.0 / 3.0)

    return radius, shell


@njit(**njit_dict_no_parallel)
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


@njit(**njit_dict_no_parallel)
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


@njit(**njit_dict_no_parallel)
def sample_energy(energy, intensity):
    """Samples energy from energy and intensity

    Parameters
    ----------
    energy :  One-dimensional Numpy Array, dtype float
        Array of energies
    intensity :  One-dimensional Numpy Array, dtype float
        Array of intensities

    Returns
    -------
    float
        Energy
    """
    z = np.random.random()

    average = (energy * intensity).sum()
    total = 0
    for e, i in zip(energy, intensity):
        total += e * i / average
        if z <= total:
            return e

    return False


@njit(**njit_dict_no_parallel)
def sample_decay_time(
    start_tau, end_tau=0.0, decay_time_min=0.0, decay_time_max=0.0
):
    """Samples the decay time from the mean lifetime
    of the isotopes (needs restructuring for more isotopes)

    Parameters
    ----------
    start_tau : float64
        Initial isotope mean lifetime
    end_tau : float64, optional
        Ending mean lifetime, by default 0 for single decays

    Returns
    -------
    float64
        Sampled decay time
    """
    decay_time = decay_time_min
    while (decay_time <= decay_time_min) or (decay_time >= decay_time_max):
        decay_time = -start_tau * np.log(np.random.random()) - end_tau * np.log(
            np.random.random()
        )
    return decay_time


class PositroniumSampler:
    def __init__(self, n_grid=10000):
        """
        Parameters
        ----------
        n_grid : int, optional
            Number of grid points for the CDF, by default 1000
        """
        self.x_grid = np.linspace(1e-4, 0.9999, n_grid)
        self.norm_pdf = self.pdf(self.x_grid) / np.trapz(self.pdf(self.x_grid), self.x_grid)
        self.cdf_grid = np.cumsum(self.norm_pdf)
        self.cdf_grid /= self.cdf_grid[-1]

    @staticmethod
    def pdf(x):
        """
        Parameters
        ----------
        x : float
            E / m_e c^2

        Returns
        -------
        float
            PDF of positronium energy from Ore amnd Powell (1949)
        """
        first_term = x * (1 - x) / (2 - x) ** 2
        second_term = 2 * (1 - x) ** 2 * np.log(1 - x) / (2 - x) ** 3
        third_term = (2 - x) / x
        fourth_term = 2 * (1 - x) * np.log(1 - x) / x**2

        return 2 * (first_term - second_term + third_term + fourth_term)
    
    def sample_energy(self, samples):
        """
        Returns
        -------
        float
            Sampled positronium energy
        """
        inverse_cdf = interp1d(self.cdf_grid, self.x_grid, bounds_error=False, fill_value="extrapolate")

        z = np.random.random(samples)
        # converted to keV
        return (inverse_cdf(z) * const.m_e.cgs.value 
            * const.c.cgs.value**2 
         * u.erg.to(u.keV))