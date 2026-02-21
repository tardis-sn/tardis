import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy import units as u
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z

from tardis.io.configuration.config_reader import Configuration
from tardis.io.configuration.config_validator import validate_dict
from tardis.io.model.csvy.data import CSVYData
from tardis.io.util import YAMLLoader
from tardis.util.base import is_valid_nuclide_or_elem, quantity_linspace

YAML_DELIMITER = "---"

logger = logging.getLogger(__name__)


def load_csvy(fname: str | Path) -> CSVYData:
    """
    Load CSVY file and return a CSVYData dataclass object.

    This function extracts velocity, density, and mass fractions from the CSVY file
    and constructs a CSVYData object. Velocity and density can be defined either
    in the YAML metadata or in the CSV data section.

    Parameters
    ----------
    fname : str or pathlib.Path
        Path to csvy file.

    Returns
    -------
    CSVYData
        Dataclass containing YAML metadata, velocity, density, and mass fractions.
    """
    fname = Path(fname)
    with fname.open() as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert line.strip() == YAML_DELIMITER, (
                    "First line of csvy file is not '---'"
                )
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
        try:
            csvy_data = pd.read_csv(fname, skiprows=yaml_end_ind + 1)
        except pd.errors.EmptyDataError as e:
            logger.debug("Could not Read CSV. Setting Dataframe to None")
            csvy_data = None

    # Extract velocity
    if "velocity" in yaml_dict and isinstance(yaml_dict["velocity"], dict):
        # Velocity defined in YAML as linspace parameters
        velocity = quantity_linspace(
            yaml_dict["velocity"]["start"],
            yaml_dict["velocity"]["stop"],
            yaml_dict["velocity"]["num"] + 1,
        ).cgs.value
    elif csvy_data is not None and "velocity" in csvy_data.columns:
        # Velocity defined in CSV data
        velocity_field_index = [
            field["name"] for field in yaml_dict["datatype"]["fields"]
        ].index("velocity")
        velocity_unit = u.Unit(
            yaml_dict["datatype"]["fields"][velocity_field_index]["unit"]
        )
        velocity = (
            (csvy_data["velocity"].values * velocity_unit).to("cm/s").value
        )
    else:
        raise ValueError("Velocity information not found in CSVY file")

    # Extract density (similar logic)
    if "density" in yaml_dict and isinstance(yaml_dict["density"], dict):
        # Density is parametric - will be calculated later, set to None for now
        density = None
    elif csvy_data is not None and "density" in csvy_data.columns:
        # Density defined in CSV data
        density_field_index = [
            field["name"] for field in yaml_dict["datatype"]["fields"]
        ].index("density")
        density_unit = u.Unit(
            yaml_dict["datatype"]["fields"][density_field_index]["unit"]
        )
        density = (
            (csvy_data["density"].values * density_unit).to("g/cm^3").value
        )
    else:
        density = None

    # Extract mass fractions if present in CSV
    mass_fractions = pd.DataFrame()
    isotope_mass_fractions = pd.DataFrame()

    if csvy_data is not None:
        mass_fraction_cols = [
            col for col in csvy_data.columns if is_valid_nuclide_or_elem(col)
        ]
        if mass_fraction_cols:
            _, mass_fractions, isotope_mass_fractions = (
                parse_csv_mass_fractions(csvy_data)
            )
            # Remove first column (velocity column index)
            if not mass_fractions.empty:
                mass_fractions = mass_fractions.iloc[:, 1:]
                mass_fractions.columns = np.arange(mass_fractions.shape[1])
            if not isotope_mass_fractions.empty:
                isotope_mass_fractions = isotope_mass_fractions.iloc[:, 1:]
                isotope_mass_fractions.columns = np.arange(
                    isotope_mass_fractions.shape[1]
                )

    # Validate and wrap yaml_dict into Configuration object
    csvy_schema_fname = (
        Path(__file__).parent.parent.parent
        / "configuration"
        / "schemas"
        / "csvy_model.yml"
    ).resolve()
    model_config = Configuration(
        validate_dict(yaml_dict, schemapath=csvy_schema_fname)
    )

    return CSVYData(
        model_config=model_config,
        velocity=velocity,
        density=density,
        mass_fractions=mass_fractions,
        isotope_mass_fractions=isotope_mass_fractions,
        raw_csv_data=csvy_data,
    )


def load_yaml_from_csvy(fpath: str | Path) -> dict:
    """
    Load only the YAML metadata from a CSVY file.

    Parameters
    ----------
    fname : str or pathlib.Path
        Path to csvy file.

    Returns
    -------
    dict
        YAML metadata dictionary from the csvy file.

    Raises
    ------
    AssertionError
        If the first line is not the YAML delimiter '---'.
    ValueError
        If the closing YAML delimiter '---' is not found.
    """
    fpath = Path(fpath)
    with fpath.open() as fh:
        yaml_lines = []
        for i, line in enumerate(fh):
            if i == 0:
                if line.strip() != YAML_DELIMITER:
                    raise ValueError(
                        f"First line of CSVY file must be '{YAML_DELIMITER}'"
                    )
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
                return yaml_dict
        raise ValueError(f"End {YAML_DELIMITER} not found")



def load_csv_from_csvy(fpath: str | Path) -> pd.DataFrame | None:
    """
    Load only the CSV data from a CSVY file.

    Parameters
    ----------
    fname : str or pathlib.Path
        Path to csvy file.

    Returns
    -------
    data : pandas.DataFrame or None
        CSV data from csvy file, or None if no CSV data is present.
    """
    fpath = Path(fpath)
    with fpath.open() as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert line.strip() == YAML_DELIMITER, (
                    "First line of csvy file is not '---'"
                )
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        try:
            data = pd.read_csv(fpath, skiprows=yaml_end_ind + 1)
        except pd.errors.EmptyDataError as e:
            logger.debug("Could not Read CSV. Setting Dataframe to None")
            data = None
    return data


def parse_csv_mass_fractions(
    csvy_data: pd.DataFrame,
) -> tuple[pd.Index, pd.DataFrame, pd.DataFrame]:
    """
    Parse the CSV data part of a CSVY model file and extract mass fractions.

    This function filters out columns that are not mass fractions and separates
    elemental and isotopic mass fractions.

    Parameters
    ----------
    csvy_data : pandas.DataFrame
        CSV data from CSVY file containing mass fraction columns.

    Returns
    -------
    index : pandas.Index
        Index of atomic numbers for elemental mass fractions.
    mass_fractions : pandas.DataFrame
        DataFrame of elemental mass fractions with atomic_number as index.
    isotope_mass_fractions : pandas.DataFrame
        DataFrame of isotopic mass fractions with MultiIndex of
        (atomic_number, mass_number).
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
