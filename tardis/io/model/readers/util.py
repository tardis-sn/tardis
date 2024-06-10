import numpy as np
import pandas as pd
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z


def read_csv_isotope_mass_fractions(
    fname, delimiter=r"\s+", skip_columns=0, skip_rows=[1]
):
    """
    A generic parser for a TARDIS composition stored as a CSV file

    The parser can read in both elemental and isotopic mass fractions. The first
    column is always expected to contain a running index, labelling the grid
    cells. The parser also allows for additional information to be stored in
    the first skip_columns columns. These will be ignored if skip_columns > 0.
    Note that the first column, containing the cell index is not taken into
    account here.

    Specific header lines can be skipped by the skip_rows keyword argument

    It is expected that the first row of the date block (after skipping the
    rows specified in skip_rows) specifies the different elements and isotopes.
    Each row after contains the composition in the corresponding grid shell.
    The first composition row describes the composition of the photosphere and
    is essentially ignored (for the default value of skip_rows).

    Example:

    Index C   O   Ni56
    0     1   1   1
    1     0.4 0.3 0.2

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    index : np.ndarray
    mass_fractions : pandas.DataFrame
    isotope_mass_fraction : pandas.MultiIndex
    """
    df = pd.read_csv(
        fname, comment="#", sep=delimiter, skiprows=skip_rows, index_col=0
    )
    df = df.transpose()

    mass_fractions = pd.DataFrame(
        columns=np.arange(df.shape[1]),
        index=pd.Index([], name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_mass_fractions = pd.DataFrame(
        columns=np.arange(df.shape[1]), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in df.index[skip_columns:]:
        if element_symbol_string in Z_DICT.values():
            z = elem_to_Z(element_symbol_string)
            mass_fractions.loc[z, :] = df.loc[element_symbol_string].tolist()
        else:
            nuc = Nuclide(element_symbol_string)
            z = nuc.Z
            mass_no = nuc.A
            isotope_mass_fractions.loc[(z, mass_no), :] = df.loc[
                element_symbol_string
            ].tolist()

    return mass_fractions.index, mass_fractions, isotope_mass_fractions
