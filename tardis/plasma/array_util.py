"""Shared helpers for continuum-related plasma data structures."""

import numpy as np
import pandas as pd
from numba import njit, prange

from tardis.transport.montecarlo import njit_dict


def get_ion_multi_index(
    multi_index_full: pd.MultiIndex, next_higher: bool = True
) -> pd.MultiIndex:
    """Return the multi-index for all ions.

    Parameters
    ----------
    multi_index_full : pd.MultiIndex
        A full multi-index with at least atomic number, ion number indices.
    next_higher : bool, optional
        Take the next higher ions, by default True

    Returns
    -------
    pd.MultiIndex
        multi-index of an atom
    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


def get_ground_state_multi_index(multi_index_full: pd.MultiIndex) -> pd.MultiIndex:
    """Return the next-ion ground-state index for a level MultiIndex."""
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1) + 1
    level_number = np.zeros_like(ion_number)
    return pd.MultiIndex.from_arrays([atomic_number, ion_number, level_number])


def cooling_rate_series2dataframe(
    cooling_rate_series: pd.Series, destination_level_idx: str
) -> pd.DataFrame:
    """Format a cooling-rate Series for transition-probability collection.

    Returns
    -------
    pd.DataFrame
        Cooling rates as a dataframe.
    """
    index = pd.MultiIndex.from_tuples(
        [("k", destination_level_idx, -1)],
        names=[
            "source_level_idx",
            "destination_level_idx",
            "transition_type",
        ],
    )
    return pd.DataFrame(cooling_rate_series.values[np.newaxis], index=index)


@njit(**njit_dict)
def numba_cumulative_trapezoid(f, x):
    """Cumulatively integrate and normalize one array using trapezoids.

    Parameters
    ----------
    f : np.ndarray
        Array to integrate
    x : np.ndarray
        Coordinate to integrate over

    Returns
    -------
    np.ndarray
        Cumulative normalized integral
    """
    integ = (np.diff(x) * (f[1:] + f[:-1]) / 2.0).cumsum()
    return integ / integ[-1]


@njit(**njit_dict)
def cumulative_integrate_array_by_blocks(f, x, block_references):
    """Compute normalized cumulative trapezoidal integrals by block.

    Parameters
    ----------
    f : np.ndarray
        Two-dimensional array of integrand values.
    x : np.ndarray
        Integration coordinates.
    block_references : np.ndarray
        Indices defining the block boundaries.

    Returns
    -------
    np.ndarray
        Normalized cumulative integrals with the same shape as ``f``.
    """
    n_rows = len(block_references) - 1
    integrated = np.zeros_like(f)
    for i in prange(f.shape[1]):
        for j in range(n_rows):
            start = block_references[j]
            stop = block_references[j + 1]
            if stop - start <= 1:
                continue
            cumulative_integral = 0.0
            for k in range(start + 1, stop):
                cumulative_integral += (
                    (x[k] - x[k - 1]) * (f[k, i] + f[k - 1, i]) / 2.0
                )
                integrated[k, i] = cumulative_integral
            for k in range(start + 1, stop):
                integrated[k, i] /= cumulative_integral
    return integrated
