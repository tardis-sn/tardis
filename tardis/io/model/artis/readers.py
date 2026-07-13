from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.model.artis.data import ArtisData


def _remove_explicit_isotopes_from_elemental_mass_fractions(
    elemental_mass_fractions: pd.DataFrame,
    isotope_mass_fractions: pd.DataFrame,
) -> pd.DataFrame:
    """Remove explicitly tracked ARTIS isotopes from elemental totals."""
    isotope_element_mass_fractions = isotope_mass_fractions.groupby(
        level="atomic_number"
    ).sum()
    residual_elemental_mass_fractions = elemental_mass_fractions.sub(
        isotope_element_mass_fractions, fill_value=0.0
    )

    roundoff_negative = (residual_elemental_mass_fractions < 0.0) & (
        residual_elemental_mass_fractions >= -np.finfo(np.float64).eps
    )
    return residual_elemental_mass_fractions.mask(roundoff_negative, 0.0)


def read_artis_density(
    fname: str | Path, legacy_return: bool = True
) -> (
    tuple[u.Quantity, u.Quantity, u.Quantity]
    | tuple[u.Quantity, u.Quantity, u.Quantity, pd.DataFrame]
):
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
    fname : str or pathlib.Path
        Path to the ARTIS density file.
    legacy_return : bool, optional
        If True, returns (time_of_model, velocity, mean_density).
        If False, returns (time_of_model, velocity, mean_density, isotope_mass_fractions).
        Default is True.

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
    fname = Path(fname)
    with fname.open() as fh:
        for i, line in enumerate(fh):
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
        dtype=dict.fromkeys(artis_model_columns, np.float64),
        names=artis_model_columns,
        sep=r"\s+",
    )
    assert len(artis_model) == no_of_shells, (
        f"Number of shells {len(artis_model)} does not match metadata {no_of_shells}"
    )
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

    isotope_mass_fractions.columns = artis_model["cell_index"].astype(np.int64)
    isotope_mass_fractions.columns.name = "cell_index"
    mean_density = mean_density[1:]
    isotope_mass_fractions = isotope_mass_fractions.iloc[:, 1:]

    if legacy_return:
        return time_of_model, velocity, mean_density
    return time_of_model, velocity, mean_density, isotope_mass_fractions


def read_artis_mass_fractions(
    fname: str | Path, normalize: bool = True
) -> pd.DataFrame:
    """
    Reads mass fractions from an ARTIS abundance file.

    The file must have shell indices in the first column, followed by
    columns for each element. Each row generally corresponds to one shell.

    Parameters
    ----------
    fname : str or pathlib.Path
        Path to the ARTIS abundance file.
    normalize : bool, optional
        If True, normalizes each row so the sum of mass fractions is 1.
        Default is True.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with 'atomic_number' as its index and 'cell_index'
        as columns, holding the mass fractions for each element and shell.
    """
    fname = Path(fname)
    mass_fractions_df = pd.read_csv(
        fname, comment="#", sep=r"\s+", header=None, index_col=0
    )

    mass_fractions_df.index.name = "cell_index"
    mass_fractions_df = mass_fractions_df.iloc[1:]

    if normalize:
        mass_fractions_df = mass_fractions_df.div(
            mass_fractions_df.sum(axis=1), axis=0
        )
    mass_fractions_df = mass_fractions_df.T
    mass_fractions_df.index.name = "atomic_number"
    mass_fractions_df.columns.name = "cell_index"

    return mass_fractions_df


def read_artis_composition(
    density_fname: str | Path, abundance_fname: str | Path
) -> tuple[pd.Index, pd.DataFrame, pd.DataFrame]:
    """
    Read ARTIS elemental and isotopic mass fractions.

    ARTIS stores total elemental mass fractions in the abundance file and
    selected radioactive isotope mass fractions in the density file. This
    reader removes those explicit isotopes from their elemental totals so that
    each nuclide is represented exactly once downstream.

    Parameters
    ----------
    density_fname : str or pathlib.Path
        Path to the ARTIS density file.
    abundance_fname : str or pathlib.Path
        Path to the ARTIS abundance file.

    Returns
    -------
    index : pandas.Index
        Atomic numbers for elemental mass fractions.
    mass_fractions : pandas.DataFrame
        Elemental mass fractions indexed by atomic number.
    isotope_mass_fractions : pandas.DataFrame
        Isotopic mass fractions indexed by atomic and mass number.
    """
    _, _, _, isotope_mass_fractions = cast(
        "tuple[u.Quantity, u.Quantity, u.Quantity, pd.DataFrame]",
        read_artis_density(density_fname, legacy_return=False),
    )
    mass_fractions = read_artis_mass_fractions(abundance_fname)
    mass_fractions = _remove_explicit_isotopes_from_elemental_mass_fractions(
        mass_fractions, isotope_mass_fractions
    )
    return mass_fractions.index, mass_fractions, isotope_mass_fractions


def read_artis_model(
    density_fname: str | Path, abundance_fname: str | Path
) -> ArtisData:
    """
    Read density and abundance files, combine them, and return a single model dataset.

    Parameters
    ----------
    density_fname : str or pathlib.Path
        Path to the ARTIS density file.
    abundance_fname : str or pathlib.Path
        Path to the ARTIS abundance file.

    Returns
    -------
    ArtisData
        Combined data with time_of_model, velocity, mean_density, and mass_fractions.
    """
    time_of_model, velocity, mean_density, isotope_mass_fractions = cast(
        "tuple[u.Quantity, u.Quantity, u.Quantity, pd.DataFrame]",
        read_artis_density(density_fname, legacy_return=False),
    )
    mass_fractions = read_artis_mass_fractions(abundance_fname)
    mass_fractions = _remove_explicit_isotopes_from_elemental_mass_fractions(
        mass_fractions, isotope_mass_fractions
    )
    mass_fractions.index = pd.MultiIndex.from_arrays(
        [mass_fractions.index, [-1] * len(mass_fractions.index)],
        names=["atomic_number", "mass_number"],
    )

    mass_fractions = pd.concat([mass_fractions, isotope_mass_fractions], axis=0)

    return ArtisData(
        time_of_model=time_of_model,
        velocity=velocity,
        mean_density=mean_density,
        mass_fractions=mass_fractions,
        isotope_mass_fractions=isotope_mass_fractions,
    )
