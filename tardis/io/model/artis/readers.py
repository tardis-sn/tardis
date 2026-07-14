from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.model.artis.data import ArtisData

ARTIS_MODEL_COLUMNS = (
    "cell_index",
    "velocities",
    "mean_densities_0",
    "ni56_mass_fraction",
    "co56_mass_fraction",
    "fe52_mass_fraction",
    "cr48_mass_fraction",
)
ARTIS_ISOTOPES = {
    "ni56_mass_fraction": (28, 56),
    "co56_mass_fraction": (27, 56),
    "fe52_mass_fraction": (26, 52),
    "cr48_mass_fraction": (24, 48),
}


def remove_explicit_isotopes_from_elemental_mass_fractions(
    elemental_mass_fractions: pd.DataFrame,
    isotope_mass_fractions: pd.DataFrame,
) -> pd.DataFrame:
    """Remove explicitly tracked ARTIS isotopes from elemental totals.

    ARTIS abundance files contain total elemental mass fractions, while the
    corresponding structure files separately identify selected isotopes. This
    function subtracts those isotope fractions from their elemental totals so
    that each nuclide is represented exactly once in the TARDIS composition.

    Parameters
    ----------
    elemental_mass_fractions : pandas.DataFrame
        Total elemental mass fractions. Rows are indexed by atomic number and
        columns correspond to model cells.
    isotope_mass_fractions : pandas.DataFrame
        Explicit isotope mass fractions. Rows use a ``(atomic_number,
        mass_number)`` MultiIndex and columns correspond to model cells.

    Returns
    -------
    pandas.DataFrame
        Residual elemental mass fractions with explicit isotope contributions
        removed. The index and columns follow ``elemental_mass_fractions``.

    Notes
    -----
    Isotope fractions are summed by atomic number before subtraction. Negative
    residuals no larger than machine precision are treated as floating-point
    roundoff and replaced with zero. Larger negative values are preserved for
    downstream composition validation.
    """
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


def _read_artis_structure(fname: str | Path) -> ArtisData:
    """Read density and isotope data from an ARTIS structure file."""
    fname = Path(fname)
    with fname.open() as fh:
        no_of_shells = int(fh.readline().strip())
        time_of_model = u.Quantity(float(fh.readline().strip()), "day").to("s")

    artis_model = pd.read_csv(
        fname,
        skiprows=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        dtype=dict.fromkeys(ARTIS_MODEL_COLUMNS, np.float64),
        names=ARTIS_MODEL_COLUMNS,
        sep=r"\s+",
    )
    assert len(artis_model) == no_of_shells, (
        f"Number of shells {len(artis_model)} does not match metadata {no_of_shells}"
    )

    isotope_mass_fractions = artis_model[list(ARTIS_ISOTOPES)].T.copy()
    isotope_mass_fractions.index = pd.MultiIndex.from_tuples(
        ARTIS_ISOTOPES.values(),
        names=["atomic_number", "mass_number"],
    )
    isotope_mass_fractions.columns = artis_model["cell_index"].astype(np.int64)
    isotope_mass_fractions.columns.name = "cell_index"

    return ArtisData(
        time_of_model=time_of_model,
        velocity=u.Quantity(artis_model["velocities"], "km/s").to("cm/s"),
        mean_density=u.Quantity(
            10 ** artis_model["mean_densities_0"], "g/cm^3"
        )[1:],
        isotope_mass_fractions=isotope_mass_fractions.iloc[:, 1:],
    )


def read_artis_density(
    fname: str | Path,
) -> tuple[u.Quantity, u.Quantity, u.Quantity]:
    """Read an ARTIS density file. Returns the time of model, velocity, and density.

    Parameters
    ----------
    fname : str or pathlib.Path
        Path to the ARTIS density file.

    Returns
    -------
    time_of_model : astropy.units.Quantity
        The time at which the model is valid, in seconds.
    velocity : astropy.units.Quantity
        The velocity array in cm/s.
    mean_density : astropy.units.Quantity
        The array of mean densities in g/cm^3, excluding the first (central) value.
    """
    artis_data = _read_artis_structure(fname)
    return (
        artis_data.time_of_model,
        artis_data.velocity,
        artis_data.mean_density,
    )


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


def _read_artis_abundances(
    abundance_fname: str | Path,
    isotope_mass_fractions: pd.DataFrame,
) -> pd.DataFrame:
    elemental_mass_fractions = read_artis_mass_fractions(abundance_fname)
    return remove_explicit_isotopes_from_elemental_mass_fractions(
        elemental_mass_fractions, isotope_mass_fractions
    )


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
    artis_data = _read_artis_structure(density_fname)
    elemental_mass_fractions = _read_artis_abundances(
        abundance_fname, artis_data.isotope_mass_fractions
    )
    return (
        elemental_mass_fractions.index,
        elemental_mass_fractions,
        artis_data.isotope_mass_fractions,
    )


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
    artis_data = _read_artis_structure(density_fname)
    elemental_mass_fractions = _read_artis_abundances(
        abundance_fname, artis_data.isotope_mass_fractions
    )
    elemental_mass_fractions.index = pd.MultiIndex.from_arrays(
        [
            elemental_mass_fractions.index,
            [-1] * len(elemental_mass_fractions.index),
        ],
        names=["atomic_number", "mass_number"],
    )
    artis_data.mass_fractions = pd.concat(
        [elemental_mass_fractions, artis_data.isotope_mass_fractions]
    )
    return artis_data
