# reading different model files

import warnings
import numpy as np
from numpy import recfromtxt, genfromtxt
import pandas as pd
from astropy import units as u
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z

import logging

# Adding logging support
logger = logging.getLogger(__name__)

from tardis.util.base import parse_quantity, is_valid_nuclide_or_elem


class ConfigurationError(Exception):
    pass


def read_density_file(filename, filetype):
    """
    read different density file formats

    Parameters
    ----------
    filename : str
        filename or path of the density file
    filetype : str
        type of the density file

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : np.ndarray
        the array containing the velocities
    unscaled_mean_densities : np.ndarray
        the array containing the densities
    """
    file_parsers = {
        "artis": read_artis_density,
        "simple_ascii": read_simple_ascii_density,
        "cmfgen_model": read_cmfgen_density,
    }

    electron_densities = None
    temperature = None
    if filetype == "cmfgen_model":
        (
            time_of_model,
            velocity,
            unscaled_mean_densities,
            electron_densities,
            temperature,
        ) = read_cmfgen_density(filename)
    else:
        (time_of_model, velocity, unscaled_mean_densities) = file_parsers[
            filetype
        ](filename)

    v_inner = velocity[:-1]
    v_outer = velocity[1:]

    invalid_volume_mask = (v_outer - v_inner) <= 0
    if invalid_volume_mask.sum() > 0:
        message = "\n".join(
            [
                f"cell {i:d}: v_inner {v_inner_i:s}, v_outer " f"{v_outer_i:s}"
                for i, v_inner_i, v_outer_i in zip(
                    np.arange(len(v_outer))[invalid_volume_mask],
                    v_inner[invalid_volume_mask],
                    v_outer[invalid_volume_mask],
                )
            ]
        )
        raise ConfigurationError(
            "Invalid volume of following cell(s):\n" f"{message:s}"
        )

    return (
        time_of_model,
        velocity,
        unscaled_mean_densities,
        electron_densities,
        temperature,
    )


def read_abundances_file(
    abundance_filename,
    abundance_filetype,
    inner_boundary_index=None,
    outer_boundary_index=None,
):
    """
    read different density file formats

    Parameters
    ----------
    abundance_filename : str
        filename or path of the density file
    abundance_filetype : str
        type of the density file
    inner_boundary_index : int
        index of the inner shell, default None
    outer_boundary_index : int
        index of the outer shell, default None
    """

    file_parsers = {
        "simple_ascii": read_simple_ascii_abundances,
        "artis": read_simple_ascii_abundances,
        "cmfgen_model": read_cmfgen_composition,
        "custom_composition": read_csv_composition,
    }

    isotope_abundance = pd.DataFrame()
    if abundance_filetype in ["cmfgen_model", "custom_composition"]:
        index, abundances, isotope_abundance = file_parsers[abundance_filetype](
            abundance_filename
        )
    else:
        index, abundances = file_parsers[abundance_filetype](abundance_filename)

    if outer_boundary_index is not None:
        outer_boundary_index_m1 = outer_boundary_index - 1
    else:
        outer_boundary_index_m1 = None
    index = index[inner_boundary_index:outer_boundary_index]
    abundances = abundances.loc[
        :, slice(inner_boundary_index, outer_boundary_index_m1)
    ]
    abundances.columns = np.arange(len(abundances.columns))
    return index, abundances, isotope_abundance


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


def read_simple_ascii_density(fname):
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
    data : pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
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


def read_artis_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5
    #index velocity [km/s] log10(density) [log10(g/cm^3)]
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
    data : pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
    """

    with open(fname) as fh:
        for i, line in enumerate(open(fname)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), "day").to("s")
            elif i == 2:
                break

    artis_model_columns = [
        "index",
        "velocities",
        "mean_densities_0",
        "ni56_fraction",
        "co56_fraction",
        "fe52_fraction",
        "cr48_fraction",
    ]
    artis_model = recfromtxt(
        fname,
        skip_header=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        unpack=True,
        dtype=[(item, np.float64) for item in artis_model_columns],
    )

    velocity = u.Quantity(artis_model["velocities"], "km/s").to("cm/s")
    mean_density = u.Quantity(10 ** artis_model["mean_densities_0"], "g/cm^3")[
        1:
    ]

    return time_of_model, velocity, mean_density


def read_cmfgen_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    The file consists of a header row and next row contains unit of the respective attributes
    Note that the first column has to contain a running index

    Example:

    index velocity densities electron_densities temperature
    - km/s g/cm^3 /cm^3 K
    0 871.66905 4.2537191e-09 2.5953807e+14 7.6395577
    1 877.44269 4.2537191e-09 2.5953807e+14 7.6395577

    Rest columns contain abundances of elements and isotopes

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : np.ndarray
    mean_density : np.ndarray
    electron_densities : np.ndarray
    temperature : np.ndarray
    """
    warnings.warn(
        "The current CMFGEN model parser is deprecated", DeprecationWarning
    )

    df = pd.read_csv(fname, comment="#", delimiter=r"\s+", skiprows=[0, 2])

    with open(fname) as fh:
        for row_index, line in enumerate(fh):
            if row_index == 0:
                time_of_model_string = line.strip().replace("t0:", "")
                time_of_model = parse_quantity(time_of_model_string)
            elif row_index == 2:
                quantities = line.split()

    velocity = u.Quantity(df["velocity"].values, quantities[1]).to("cm/s")
    temperature = u.Quantity(df["temperature"].values, quantities[2])[1:]
    mean_density = u.Quantity(df["densities"].values, quantities[3])[1:]
    electron_densities = u.Quantity(
        df["electron_densities"].values, quantities[4]
    )[1:]

    return (
        time_of_model,
        velocity,
        mean_density,
        electron_densities,
        temperature,
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


def read_cmfgen_composition(fname, delimiter=r"\s+"):
    """Read composition from a CMFGEN model file

    The CMFGEN file format contains information about the ejecta state in the
    first four columns and the following ones contain elemental and isotopic
    abundances.

    WARNING : deprecated

    fname : str
        filename of the csv file
    """

    warnings.warn(
        "The current CMFGEN model parser is deprecated", DeprecationWarning
    )

    return read_csv_isotope_abundances(
        fname, delimiter=delimiter, skip_columns=4, skip_rows=[0, 2, 3]
    )


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


def read_csv_isotope_abundances(
    fname, delimiter=r"\s+", skip_columns=0, skip_rows=[1]
):
    """
    A generic parser for a TARDIS composition stored as a CSV file

    The parser can read in both elemental and isotopic abundances. The first
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
    abundances : pandas.DataFrame
    isotope_abundance : pandas.MultiIndex
    """

    df = pd.read_csv(
        fname, comment="#", sep=delimiter, skiprows=skip_rows, index_col=0
    )
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

    for element_symbol_string in df.index[skip_columns:]:
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
