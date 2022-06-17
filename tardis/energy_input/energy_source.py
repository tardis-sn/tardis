import numpy as np
import pandas as pd
from numba import njit
import radioactivedecay as rd

from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.util.base import (
    atomic_number2element_symbol,
)
from tardis.energy_input.util import (
    convert_half_life_to_astropy_units,
)


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
    float
        Number of decay products per decay
    """
    intensity_1 = nuclear_data.query("type==" + source_1)["intensity"].values
    intensity_2 = nuclear_data.query("type==" + source_2)["intensity"].values
    total_intensity = np.sum(intensity_1) + np.sum(intensity_2)
    # Factor of 100 is because intensities are in percent
    scale_factor = total_intensity / 100
    return (
        np.sum(intensity_1) / total_intensity,
        np.sum(intensity_2) / total_intensity,
        scale_factor,
    )


def decay_chain_energy(
    time_start,
    time_end,
    number_start_isotope,
    start_isotope_tau,
    end_isotope_tau,
    start_isotope_energy,
    end_isotope_energy,
    start_isotope_name,
    end_isotope_name,
):
    """Calculate the energy from an arbitrary decay chain

    Parameters
    ----------
    time_start : float
        Start time in days
    time_end : float
        End time in days
    number_start_isotope : int
        Number of initial isotope atoms at time_start
    start_isotope_tau, end_isotope_tau : float64
        Mean half-life for each isotope
    start_isotope_energy : float64
        Initial isotope total energy
    end_isotope_energy : float64
        Final isotope total energy

    Returns
    -------
    DataFrame
        Total energy from decay for each isotope
    """
    total_start_isotope = -start_isotope_tau * (
        np.exp(-time_end / start_isotope_tau)
        - np.exp(-time_start / start_isotope_tau)
    )
    total_end_isotope = -end_isotope_tau * (
        np.exp(-time_end / end_isotope_tau)
        - np.exp(-time_start / end_isotope_tau)
    )

    total_energy = pd.DataFrame()

    total_energy[start_isotope_name] = number_start_isotope * (
        (start_isotope_energy / start_isotope_tau) * total_start_isotope
    )

    total_energy[end_isotope_name] = number_start_isotope * (
        (end_isotope_energy / (start_isotope_tau - end_isotope_tau))
        * (total_start_isotope - total_end_isotope)
    )

    return total_energy


def get_all_isotopes(abundances):
    """Get the possible isotopes present over time
    for a given starting abundance

    Parameters
    ----------
    abundances : DataFrame
        Current isotope abundances

    Returns
    -------
    list
        List of isotope names
    """
    progenitors = [
        f"{rd.utils.Z_DICT[i[0]]}-{i[1]}" for i in abundances.T.columns
    ]

    isotopes = set(progenitors)
    check = True

    while check == True:
        progeny = set(isotopes)

        for i in isotopes:
            for p in rd.Nuclide(i).progeny():
                if (
                    p != "SF"
                    and rd.Nuclide(p).half_life("readable") != "stable"
                ):
                    progeny.add(p)

        if progeny == isotopes:
            check = False
        else:
            isotopes |= progeny

    isotopes = [i.replace("-", "") for i in isotopes]
    return isotopes


def get_progeny_for_isotopes(abundances):
    """Produces a nested list of all direct radioactive progeny
    of the input abundance.

    Parameters
    ----------
    abundances : DataFrame
        Current isotope abundances

    Returns
    -------
    List
        List of progeny lists indexed matching the abundance columns
    """
    progenitors = [
        f"{rd.utils.Z_DICT[i[0]]}-{i[1]}" for i in abundances.T.columns
    ]

    progeny = []

    for isotope in progenitors:
        isotope_progeny = []
        for p in rd.Nuclide(isotope).progeny():
            if p != "SF" and rd.Nuclide(p).half_life("readable") != "stable":
                isotope_progeny.append(p.replace("-", ""))

        progeny.append(isotope_progeny)

    return progeny


def get_tau(meta, isotope_string):
    """Calculate the mean lifetime of an isotope

    Parameters
    ----------
    meta : DataFrame
        Isotope metadata
    isotope_string : str
        Isotope of interest

    Returns
    -------
    float
        Mean lifetime of isotope
    """
    isotope_meta = meta.loc[isotope_string]
    half_life = isotope_meta.loc[isotope_meta["key"] == "Parent T1/2 value"][
        "value"
    ].values[0]
    half_life = convert_half_life_to_astropy_units(half_life)
    return half_life / np.log(2)


def get_isotope_string(atom_number, atom_mass):
    """Get the isotope string in the format e.g. Ni56

    Parameters
    ----------
    atom_number : int
        Atomic number
    atom_mass : int
        Atomic mass

    Returns
    -------
    str
        Isotope string in the format e.g. Ni56
    """
    return atomic_number2element_symbol(atom_number) + str(atom_mass)


def read_artis_lines(isotope, path_to_data):
    """Reads lines of ARTIS format

    Parameters
    ----------
    isotope : string
        Isotope to read e.g. Ni56

    Returns
    -------
    pd.DataFrame
        Energies and intensities of the isotope lines
    """
    return pd.read_csv(
        path_to_data + isotope + ".txt",
        names=["energy", "intensity"],
        sep="  ",
        index_col=False,
    )
