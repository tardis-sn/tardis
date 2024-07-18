import logging
import numpy as np
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z
import yaml
import pandas as pd
from tardis.io.util import YAMLLoader
from tardis.util.base import is_valid_nuclide_or_elem

YAML_DELIMITER = "---"

logger = logging.getLogger(__name__)


def load_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    yaml_dict : dictionary
                YAML part of the csvy file
    data : pandas.dataframe
            csv data from csvy file
    """
    with open(fname) as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert (
                    line.strip() == YAML_DELIMITER
                ), "First line of csvy file is not '---'"
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
        try:
            data = pd.read_csv(fname, skiprows=yaml_end_ind + 1)
        except pd.errors.EmptyDataError as e:
            logger.debug("Could not Read CSV. Setting Dataframe to None")
            data = None

    return yaml_dict, data


def load_yaml_from_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    yaml_dict : dictionary
                YAML part of the csvy file
    """
    with open(fname) as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert (
                    line.strip() == YAML_DELIMITER
                ), "First line of csvy file is not '---'"
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
    return yaml_dict


def load_csv_from_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    data : pandas.dataframe
           csv data from csvy file
    """
    yaml_dict, data = load_csvy(fname)
    return data


def parse_csv_mass_fractions(csvy_data):
    """
    A parser for the csv data part of a csvy model file. This function filters out columns that are not mass fractions.

    Parameters
    ----------
    csvy_data : pandas.DataFrame

    Returns
    -------
    index : np.ndarray
    mass_fractions : pandas.DataFrame
    isotope_mass_fraction : pandas.MultiIndex
    """

    mass_fraction_col_names = [
        name for name in csvy_data.columns if is_valid_nuclide_or_elem(name)
    ]
    df = csvy_data.loc[:, mass_fraction_col_names]

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

    for element_symbol_string in df.index[0:]:
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
