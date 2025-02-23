from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u


@dataclass
class ArtisModelData:
    time_of_model: u.Quantity
    velocity: np.ndarray
    mean_density: u.Quantity
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def to_geometry(self):
        """
        Construct a HomologousRadial1DGeometry object from this ArtisModelData.

        We create v_inner and v_outer by treating the velocity array as boundary
        points for the shells. The time_of_model is used as the time_explosion.
        """
        geometry = HomologousRadial1DGeometry(
            v_inner=self.velocity[:-1],
            v_outer=self.velocity[1:],
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=self.time_of_model,
        )
        return geometry


def read_artis_density(fname, legacy_return=True):
    """
    Read an ARTIS density file.

    The file is expected to have:
      • First line: number of shells (ignored beyond reading the integer).
      • Second line: time of the model in days, which is then converted to seconds.
      • Remaining lines: columns containing cell_id, velocity (km/s), log10(density),
        and mass fractions of Ni56, Co56, Fe52, and Cr48.

    Notes
    -----
    The first density (i.e., the center of the model) is not used in the returned arrays.

    Parameters
    ----------
    fname : str
        Path to the ARTIS density file.
    legacy_return : bool, optional (default: True)
        If True, returns (time_of_model, velocity, mean_density).
        If False, returns (time_of_model, velocity, mean_density, isotope_mass_fractions).

    Returns
    -------
    time_of_model : astropy.units.Quantity
        The time at which the model is valid, in seconds.
    velocity : astropy.units.Quantity
        The velocity array in cm/s.
    mean_density : astropy.units.Quantity
        The array of mean densities in g/cm^3, excluding the first (central) value.
    isotope_mass_fractions : pandas.DataFrame, optional
        Mass fractions for Ni56, Co56, Fe52, and Cr48 if legacy_return=False.
        The DataFrame has a MultiIndex of (Z, A) for rows and the cell_id for columns.
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
        "cell_index",
        "velocities",
        "mean_densities_0",
        "ni56_mass_fraction",
        "co56_mass_fraction",
        "fe52_mass_fraction",
        "cr48_mass_fraction",
    ]

    artis_model = pd.read_csv(
        fname,
        skiprows=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        dtype={item: np.float64 for item in artis_model_columns},
        names=artis_model_columns,
        sep=r"\s+",
    )
    assert (
        len(artis_model) == no_of_shells
    ), "Number of shells {len(artis_model)} does not match metadate {no_of_shells}"
    velocity = u.Quantity(artis_model["velocities"], "km/s").to("cm/s")
    mean_density = u.Quantity(10 ** artis_model["mean_densities_0"], "g/cm^3")

    isotope_mass_fractions = pd.DataFrame(
        artis_model[
            [
                "ni56_mass_fraction",
                "co56_mass_fraction",
                "fe52_mass_fraction",
                "cr48_mass_fraction",
            ]
        ].T
    )
    isotope_mass_fractions.index = pd.MultiIndex.from_tuples(
        [(28, 56), (27, 56), (26, 52), (24, 48)],
        names=["atomic_number", "mass_number"],
    )

    isotope_mass_fractions.columns = artis_model["cell_index"]
    isotope_mass_fractions.columns.name = "cell_index"

    if legacy_return:
        return time_of_model, velocity, mean_density[1:]
    else:
        return time_of_model, velocity, mean_density, isotope_mass_fractions


def read_artis_mass_fractions(fname, normalize=True):
    """
    Reads mass fractions from an ARTIS abundance file named e.g. 'artis_abundances.dat'.
    Each row typically corresponds to one shell. The first column is shell index,
    followed by columns for each element.

    Parameters
    ----------
    fname : str
        filename or path to ARTIS abundances

    Returns
    -------
    mass_fractions : pandas.DataFrame
        DataFrame containing mass fractions per shell (rows).
        The first column is assumed to be shell index.
    """
    # Skip comment lines or set them to '#' if needed
    # Here we assume no header row, so we'll assign column names ourselves
    mass_fractions_df = pd.read_csv(
        fname, comment="#", sep="\s+", header=None, index_col=0
    )

    mass_fractions_df.index.name = "cell_index"

    if normalize:
        mass_fractions_df = mass_fractions_df.div(mass_fractions_df.sum(axis=1), axis=0)
    mass_fractions_df = mass_fractions_df.T
    mass_fractions_df.index.name = "atomic_number"
    mass_fractions_df.columns.name = "cell_index"

    return mass_fractions_df


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
    time_of_model, velocity, mean_density, isotope_mass_fractions = read_artis_density(
        density_fname, legacy_return=False
    )
    mass_fractions = read_artis_mass_fractions(abundance_fname)
    mass_fractions.index = pd.MultiIndex.from_arrays(
        [mass_fractions.index, [-1] * len(mass_fractions.index)],
        names=["atomic_number", "mass_number"],
    )
    isotope_summed = isotope_mass_fractions.groupby(level="atomic_number").sum()
    mass_fractions = mass_fractions.sub(
        isotope_summed, level="atomic_number", fill_value=0
    )

    mass_fractions = pd.concat([mass_fractions, isotope_mass_fractions], axis=0)

    # Build the dataclass with the combined info:
    return ArtisModelData(
        time_of_model=time_of_model,
        velocity=velocity,  # or keep as Quantity if you prefer
        mean_density=mean_density,
        mass_fractions=mass_fractions,
    )


from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
