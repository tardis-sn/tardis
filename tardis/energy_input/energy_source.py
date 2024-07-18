import numpy as np
import pandas as pd
import radioactivedecay as rd

from tardis.util.base import (
    atomic_number2element_symbol,
)
from tardis.energy_input.util import (
    convert_half_life_to_astropy_units,
    ELECTRON_MASS_ENERGY_KEV,
)


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
    intensity = nuclear_data[nuclear_data.Radiation == source]["Rad Intensity"]
    energy = nuclear_data[nuclear_data.Radiation == source]["Rad Energy"]

    intensity /= 100  # = [i / 100 if i > 1 else i for i in intensity]

    return energy, intensity


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

    isotopes = [i for i in isotopes]
    return isotopes


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


def get_nuclear_lines_database(
    path,
):
    """Load the nuclear decay line data set

    Parameters
    ----------
    path : str
        Path to the data set HDF file

    Returns
    -------
    pandas DataFrame
        The decay radiation lines
    """
    decay_radiation_db = pd.read_hdf(path, "decay_radiation_data")
    return decay_radiation_db


def positronium_continuum():
    """Produces a continuum of positronium decay energy
    using the function defined by Ore and Powell 1949
    and adapted by Leung 2022 to be in terms of electron
    rest mass energy

    Returns
    -------
    energy
        An array of photon energies in keV
    intensity
        An array of intensities between 0 and 1
    """

    energy = np.linspace(1, ELECTRON_MASS_ENERGY_KEV, num=100, endpoint=False)

    x = energy / ELECTRON_MASS_ENERGY_KEV

    one_minus_x = 1 - x

    term_1 = (x * one_minus_x) / (2 - x) ** 2
    term_2 = (2 * one_minus_x**2) / (2 - x) ** 3 * np.log(one_minus_x)
    term_3 = (2 - x) / x
    term_4 = (2 * one_minus_x) / x**2 * np.log(one_minus_x)

    intensity = 2 * (term_1 - term_2 + term_3 + term_4)

    return energy, intensity / np.max(intensity)
