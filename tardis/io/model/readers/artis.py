import numpy as np
import pandas as pd
from astropy import units as u
from dataclasses import dataclass, field


@dataclass
class ArtisModelData:
    time_of_model: u.Quantity
    velocity: np.ndarray
    mean_density: u.Quantity
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)


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

    artis_model = pd.read_csv(
        fname,
        skiprows=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        dtype={item: np.float64 for item in artis_model_columns},
        names=artis_model_columns,
        # The argument `delim_whitespace` was changed to `sep`
        #   because the first one is deprecated since version 2.2.0.
        #   The regular expression means: the separation is one or
        #   more spaces together (simple space, tabs, new lines).
        sep=r"\s+",
    ).to_records(index=False)

    velocity = u.Quantity(artis_model["velocities"], "km/s").to("cm/s")
    mean_density = u.Quantity(10 ** artis_model["mean_densities_0"], "g/cm^3")[1:]

    return time_of_model, velocity, mean_density


def read_artis_mass_fractions(fname):
    """
    Reads mass fractions from an artis abundance file named e.g. 'artis_abundances.dat'.
    Each row typically corresponds to one shell. The first column might be shell index,
    followed by columns for each element.

    Parameters
    ----------
    fname : str
        filename or path to artis abundances

    Returns
    -------
    mass_fractions : pandas.DataFrame
        DataFrame containing mass fractions per shell (rows).
        The first column is assumed to be shell index.
    """
    # Skip comment lines or set them to '#' if needed
    # Here we assume no header row, so we'll assign column names ourselves
    df = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None)

    # rename columns so that col0 is 'shell_index'
    col_names = [f"col_{i}" for i in range(df.shape[1])]
    df.columns = col_names
    df = df.rename(columns={"col_0": "shell_index"})

    # Convert everything except 'shell_index' into mass fractions
    # In real usage, you might want to give meaningful element names
    # or read them from a separate file. Here we build placeholders.
    ncols = df.shape[1] - 1
    element_names = [f"X_{i}" for i in range(ncols)]
    rename_map = {f"col_{i+1}": element_names[i] for i in range(ncols)}
    df = df.rename(columns=rename_map)

    # set 'shell_index' as index for clarity
    df = df.set_index("shell_index")

    return df


def read_artis_model(density_fname, abundance_fname):
    """
    Reads both the density and abundance files and returns combined data.

    Parameters
    ----------
    density_fname : str
        Path to the ARTIS density file
    abundance_fname : str
        Path to the ARTIS abundance file

    Returns
    -------
    ArtisModelData
        Contains time_of_model, velocity, mean_density, and mass_fractions
    """
    time_of_model, velocity, mean_density = read_artis_density(density_fname)
    mass_fractions = read_artis_mass_fractions(abundance_fname)

    # Build the dataclass with the combined info:
    return ArtisModelData(
        time_of_model=time_of_model,
        velocity=velocity.value,  # or keep as Quantity if you prefer
        mean_density=mean_density,
        mass_fractions=mass_fractions,
    )
