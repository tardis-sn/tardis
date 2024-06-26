from tardis.io.model.readers.util import read_csv_isotope_mass_fractions
from tardis.util.base import parse_quantity


import pandas as pd
from astropy import units as u


import warnings


def read_cmfgen_density(fname: str):
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

    return read_csv_isotope_mass_fractions(
        fname, delimiter=delimiter, skip_columns=4, skip_rows=[0, 2, 3]
    )
