from tardis.util.base import is_valid_nuclide_or_elem


import numpy as np
import pandas as pd
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z


def parse_csv_abundances(csvy_data):
    """
    A parser for the csv data part of a csvy model file. This function filters out columns that are not abundances.

    Parameters
    ----------
    csvy_data : pandas.DataFrame

    Returns
    -------
    index : np.ndarray
    abundances : pandas.DataFrame
    isotope_abundance : pandas.MultiIndex
    """

    abundance_col_names = [
        name for name in csvy_data.columns if is_valid_nuclide_or_elem(name)
    ]
    df = csvy_data.loc[:, abundance_col_names]

    df = df.transpose()

    abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]),
        index=pd.Index([], name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in df.index[0:]:
        if element_symbol_string in Z_DICT.values():
            z = elem_to_Z(element_symbol_string)
            abundance.loc[z, :] = df.loc[element_symbol_string].tolist()
        else:
            nuc = Nuclide(element_symbol_string)
            z = nuc.Z
            mass_no = nuc.A
            isotope_abundance.loc[(z, mass_no), :] = df.loc[
                element_symbol_string
            ].tolist()

    return abundance.index, abundance, isotope_abundance