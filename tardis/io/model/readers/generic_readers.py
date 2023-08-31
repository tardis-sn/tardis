from astropy import units as u
from typing import Any, Tuple
from numpy import recfromtxt
import pandas as pd
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z
from pathlib import Path

import numpy as np
from tardis.io.model.readers.util import read_csv_isotope_abundances

from tardis.util.base import parse_quantity


class ConfigurationError(Exception):
    pass


def read_simple_ascii_density(
    fname: Any,
) -> Tuple[u.Quantity, u.Quantity, u.Quantity]:
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5 s
    #index velocity [km/s] density [g/cm^3]
    0 1.1e4 1.6e8
    1 1.2e4 1.7e8

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : u.Quantity
        velocity
    mean_density: u.Quantity
        mean density
    """

    with open(fname) as fh:
        time_of_model_string = fh.readline().strip()
        time_of_model = parse_quantity(time_of_model_string)

    data = recfromtxt(
        fname,
        skip_header=1,
        names=("index", "velocity", "density"),
        dtype=(int, float, float),
    )
    velocity = (data["velocity"] * u.km / u.s).to("cm/s")
    mean_density = (data["density"] * u.Unit("g/cm^3"))[1:]

    return time_of_model, velocity, mean_density


def read_csv_composition(fname, delimiter=r"\s+"):
    """Read composition from a simple CSV file

    The CSV file can contain specific isotopes or elemental abundances in the
    different columns. The first row must contain the header in which the
    contents of each column is specified by the elemental symbol (for elemental
    abundances) or by the symbol plus mass number (for isotopic abundances).

    Example: C O Fe Ni56 Co

    The i-th row specifies the composition in the i-th shell

    fname : str
        filename of the csv file
    """

    return read_csv_isotope_abundances(
        fname, delimiter=delimiter, skip_columns=0, skip_rows=[1]
    )


def read_simple_ascii_abundances(fname):
    """
    Reading an abundance file of the following structure (example; lines starting with hash will be ignored):
    The first line of abundances describe the abundances in the center of the model and are not used.
    #index element1, element2, ..., element30
    0 0.4 0.3, .. 0.2

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    index : np.ndarray
        containing the indices
    abundances : pandas.DataFrame
        data frame containing index, element1 - element30 and columns according to the shells
    """
    data = np.loadtxt(fname)

    index = data[1:, 0].astype(int)
    abundances = pd.DataFrame(
        data[1:, 1:].transpose(), index=np.arange(1, data.shape[1])
    )

    return index, abundances


def read_uniform_abundances(abundances_section, no_of_shells):
    """
    Parameters
    ----------
    abundances_section : config.model.abundances
    no_of_shells : int

    Returns
    -------
    abundance : pandas.DataFrame
    isotope_abundance : pandas.DataFrame
    """
    abundance = pd.DataFrame(
        columns=np.arange(no_of_shells),
        index=pd.Index(np.arange(1, 120), name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_abundance = pd.DataFrame(
        columns=np.arange(no_of_shells), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in abundances_section:
        if element_symbol_string == "type":
            continue
        try:
            if element_symbol_string in Z_DICT.values():
                z = elem_to_Z(element_symbol_string)
                abundance.loc[z] = float(
                    abundances_section[element_symbol_string]
                )
            else:
                nuc = Nuclide(element_symbol_string)
                mass_no = nuc.A
                z = nuc.Z
                isotope_abundance.loc[(z, mass_no), :] = float(
                    abundances_section[element_symbol_string]
                )

        except RuntimeError as err:
            raise RuntimeError(
                f"Abundances are not defined properly in config file : {err.args}"
            )

    return abundance, isotope_abundance
